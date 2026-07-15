use std::collections::{HashMap, HashSet};
use serde::{Deserialize, Serialize};

/// Information about a single codon
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CodonInfo {
    pub aa_code: char,          // single-letter amino acid code, '*' for stop
    pub amino_acid: String,     // three-letter abbreviation
    pub codon: String,          // the 3-letter codon
    pub subfam: String,         // subfamily identifier (e.g., "Leu_CU")
}

/// A genetic code table mapping codons to amino acid info
#[derive(Debug, Clone)]
pub struct CodonTable {
    pub gcid: String,
    pub name: String,
    pub codon_map: HashMap<String, CodonInfo>,
    pub aa_to_codons: HashMap<char, Vec<String>>,
    pub start_codons: HashSet<String>,
    pub stop_codons: HashSet<String>,
    /// All 64 codons in canonical order (TTT, TTC, TTA, ...)
    pub all_codons: Vec<String>,
}

impl CodonTable {
    /// Return the standard genetic code (NCBI code 1)
    pub fn standard() -> Self {
        Self::from_ncbi_id("1").unwrap()
    }

    /// Build a codon table from an NCBI genetic code ID (1-33)
    pub fn from_ncbi_id(gcid: &str) -> Result<Self, String> {
        let (name, aa_map) = get_ncbi_code(gcid)?;
        Self::build(gcid, &name, &aa_map)
    }

    /// Build a codon table from a custom amino-acid-to-codon mapping.
    /// `aa2codon` is a list of (three_letter_aa, codon) pairs.
    pub fn from_custom(aa2codon: &[(String, String)]) -> Result<Self, String> {
        // Infer single-letter codes
        let mut custom_map: HashMap<String, char> = HashMap::new();
        for (aa3, codon) in aa2codon {
            let aa1 = aa3_to_aa1(aa3)?;
            custom_map.insert(codon.to_uppercase(), aa1);
        }
        let name = "Custom".to_string();
        Self::build("custom", &name, &custom_map)
    }

    fn build(gcid: &str, name: &str, codon_aa: &HashMap<String, char>) -> Result<Self, String> {
        let bases = ['T', 'C', 'A', 'G'];
        let mut all_codons = Vec::with_capacity(64);
        for b1 in &bases {
            for b2 in &bases {
                for b3 in &bases {
                    all_codons.push(format!("{}{}{}", b1, b2, b3));
                }
            }
        }

        let mut codon_map: HashMap<String, CodonInfo> = HashMap::new();
        let mut aa_to_codons: HashMap<char, Vec<String>> = HashMap::new();
        let mut start_codons = HashSet::new();
        let mut stop_codons = HashSet::new();

        // Process all 64 codons
        for codon in &all_codons {
            let aa1 = *codon_aa.get(codon).unwrap_or(&'X');
            let aa3 = aa1_to_aa3(aa1);

            // Determine subfamily: first two bases + amino acid
            let prefix = &codon[..2];
            let subfam = format!("{}_{}", aa3, prefix);

            let info = CodonInfo {
                aa_code: aa1,
                amino_acid: aa3.to_string(),
                codon: codon.clone(),
                subfam,
            };

            if aa1 == '*' {
                stop_codons.insert(codon.clone());
            }

            codon_map.insert(codon.clone(), info);
            aa_to_codons.entry(aa1).or_default().push(codon.clone());
        }

        // Always set standard start codons (ATG/Methionine)
        // Note: some alternative codes may have additional start codons
        if let Some(aa) = codon_aa.get("ATG") {
            if *aa == 'M' {
                start_codons.insert("ATG".to_string());
            }
        }
        // Some codes use additional start codons
        for codon in &["TTG", "CTG", "ATA", "ATT", "ATC", "GTG"] {
            if let Some(aa) = codon_aa.get(*codon) {
                if *aa == 'M' {
                    start_codons.insert(codon.to_string());
                }
            }
        }

        Ok(CodonTable {
            gcid: gcid.to_string(),
            name: name.to_string(),
            codon_map,
            aa_to_codons,
            start_codons,
            stop_codons,
            all_codons,
        })
    }

    /// Get synonymous codons for an amino acid (sorted)
    pub fn synonymous_codons(&self, aa: char) -> Vec<String> {
        let mut codons = self.aa_to_codons.get(&aa).cloned().unwrap_or_default();
        codons.sort();
        codons
    }

    /// Check if a codon is a stop codon
    pub fn is_stop(&self, codon: &str) -> bool {
        self.stop_codons.contains(codon)
    }

    /// Check if a codon is a start codon
    pub fn is_start(&self, codon: &str) -> bool {
        self.start_codons.contains(codon)
    }

    /// Get the amino acid for a codon
    pub fn get_aa(&self, codon: &str) -> Option<char> {
        self.codon_map.get(codon).map(|ci| ci.aa_code)
    }

    /// Get the subfamily for a codon
    pub fn get_subfam(&self, codon: &str) -> Option<&str> {
        self.codon_map.get(codon).map(|ci| ci.subfam.as_str())
    }

    /// List available NCBI genetic codes
    pub fn list_available() -> Vec<(String, String)> {
        NCBI_CODE_NAMES
            .iter()
            .map(|(id, name)| (id.to_string(), name.to_string()))
            .collect()
    }

    /// Group codons by subfamily
    pub fn subfamily_groups(&self) -> HashMap<String, Vec<String>> {
        let mut groups: HashMap<String, Vec<String>> = HashMap::new();
        for codon in &self.all_codons {
            if let Some(info) = self.codon_map.get(codon) {
                if info.aa_code != '*' {
                    groups
                        .entry(info.subfam.clone())
                        .or_default()
                        .push(codon.clone());
                }
            }
        }
        // Sort codons within each group
        for codons in groups.values_mut() {
            codons.sort();
        }
        groups
    }
}

/// Three-letter to one-letter amino acid code
pub fn aa3_to_aa1(aa3: &str) -> Result<char, String> {
    match aa3.to_uppercase().as_str() {
        "ALA" => Ok('A'),
        "ARG" => Ok('R'),
        "ASN" => Ok('N'),
        "ASP" => Ok('D'),
        "CYS" => Ok('C'),
        "GLN" => Ok('Q'),
        "GLU" => Ok('E'),
        "GLY" => Ok('G'),
        "HIS" => Ok('H'),
        "ILE" => Ok('I'),
        "LEU" => Ok('L'),
        "LYS" => Ok('K'),
        "MET" => Ok('M'),
        "PHE" => Ok('F'),
        "PRO" => Ok('P'),
        "SER" => Ok('S'),
        "THR" => Ok('T'),
        "TRP" => Ok('W'),
        "TYR" => Ok('Y'),
        "VAL" => Ok('V'),
        "SEC" => Ok('U'), // Selenocysteine
        "PYL" => Ok('O'), // Pyrrolysine
        "STP" | "TER" | "STOP" => Ok('*'),
        "ASX" => Ok('B'), // Asn or Asp
        "GLX" => Ok('Z'), // Gln or Glu
        "XLE" => Ok('J'), // Leu or Ile
        "XAA" | "UNK" => Ok('X'), // Unknown
        _ => Err(format!("Unknown amino acid: {aa3}")),
    }
}

/// One-letter to three-letter amino acid code
pub fn aa1_to_aa3(aa1: char) -> &'static str {
    match aa1 {
        'A' => "Ala",
        'R' => "Arg",
        'N' => "Asn",
        'D' => "Asp",
        'C' => "Cys",
        'Q' => "Gln",
        'E' => "Glu",
        'G' => "Gly",
        'H' => "His",
        'I' => "Ile",
        'L' => "Leu",
        'K' => "Lys",
        'M' => "Met",
        'F' => "Phe",
        'P' => "Pro",
        'S' => "Ser",
        'T' => "Thr",
        'W' => "Trp",
        'Y' => "Tyr",
        'V' => "Val",
        'U' => "Sec",
        'O' => "Pyl",
        '*' => "Stp",
        'B' => "Asx",
        'Z' => "Glx",
        'J' => "Xle",
        'X' => "Unk",
        _ => "Unk",
    }
}

// ================================================================
// NCBI Genetic Code Definitions
// ================================================================

/// All NCBI genetic code names
const NCBI_CODE_NAMES: &[(&str, &str)] = &[
    ("1", "Standard"),
    ("2", "Vertebrate Mitochondrial"),
    ("3", "Yeast Mitochondrial"),
    ("4", "Mold/Protozoan/Coelenterate Mitochondrial & Mycoplasma/Spiroplasma"),
    ("5", "Invertebrate Mitochondrial"),
    ("6", "Ciliate/Dasycladacean/Hexamita Nuclear"),
    ("9", "Echinoderm/Flatworm Mitochondrial"),
    ("10", "Euplotid Nuclear"),
    ("11", "Bacterial/Archaeal/Plant Plastid"),
    ("12", "Alternative Yeast Nuclear"),
    ("13", "Ascidian Mitochondrial"),
    ("14", "Alternative Flatworm Mitochondrial"),
    ("15", "Blepharisma Nuclear"),
    ("16", "Chlorophycean Mitochondrial"),
    ("21", "Trematode Mitochondrial"),
    ("22", "Scenedesmus obliquus Mitochondrial"),
    ("23", "Thraustochytrium Mitochondrial"),
    ("24", "Pterobranchia Mitochondrial"),
    ("25", "Candidate Division SR1/Gracilibacteria"),
    ("26", "Pachysolen tannophilus Nuclear"),
    ("27", "Karyorelict Nuclear"),
    ("28", "Condylostoma Nuclear"),
    ("29", "Mesodinium Nuclear"),
    ("30", "Peritrich Nuclear"),
    ("31", "Blastocrithidia Nuclear"),
    ("33", "Cephalodiscidae Mitochondrial"),
];

/// Get the NCBI genetic code by ID.
/// Returns (name, codon->aa1 mapping).
fn get_ncbi_code(gcid: &str) -> Result<(String, HashMap<String, char>), String> {
    let code_id: u32 = gcid
        .parse()
        .map_err(|_| format!("Invalid genetic code ID: {gcid}"))?;

    // Start with the standard code
    let mut mapping = STANDARD_CODE.clone();

    let name = match code_id {
        1 => "Standard",
        2 => {
            // Vertebrate Mitochondrial
            mapping.insert("AGA".into(), '*');
            mapping.insert("AGG".into(), '*');
            mapping.insert("ATA".into(), 'M');
            mapping.insert("TGA".into(), 'W');
            "Vertebrate Mitochondrial"
        }
        3 => {
            // Yeast Mitochondrial
            mapping.insert("ATA".into(), 'M');
            mapping.insert("CTT".into(), 'T');
            mapping.insert("CTC".into(), 'T');
            mapping.insert("CTA".into(), 'T');
            mapping.insert("CTG".into(), 'T');
            mapping.insert("TGA".into(), 'W');
            mapping.insert("CGA".into(), '*'); // absent, treat as stop
            mapping.insert("CGC".into(), '*'); // absent
            "Yeast Mitochondrial"
        }
        4 => {
            // Mold/Protozoan/Coelenterate Mitochondrial & Mycoplasma/Spiroplasma
            mapping.insert("TGA".into(), 'W');
            "Mold/Protozoan/Coelenterate Mitochondrial & Mycoplasma/Spiroplasma"
        }
        5 => {
            // Invertebrate Mitochondrial
            mapping.insert("AGA".into(), 'S');
            mapping.insert("AGG".into(), 'S');
            mapping.insert("ATA".into(), 'M');
            mapping.insert("TGA".into(), 'W');
            "Invertebrate Mitochondrial"
        }
        6 => {
            // Ciliate/Dasycladacean/Hexamita Nuclear
            mapping.insert("TAA".into(), 'Q');
            mapping.insert("TAG".into(), 'Q');
            "Ciliate/Dasycladacean/Hexamita Nuclear"
        }
        9 => {
            // Echinoderm/Flatworm Mitochondrial
            mapping.insert("AAA".into(), 'N');
            mapping.insert("AGA".into(), 'S');
            mapping.insert("AGG".into(), 'S');
            mapping.insert("TGA".into(), 'W');
            "Echinoderm/Flatworm Mitochondrial"
        }
        10 => {
            // Euplotid Nuclear
            mapping.insert("TGA".into(), 'C');
            "Euplotid Nuclear"
        }
        11 => {
            // Bacterial/Archaeal/Plant Plastid (same as standard)
            "Bacterial/Archaeal/Plant Plastid"
        }
        12 => {
            // Alternative Yeast Nuclear
            mapping.insert("CTG".into(), 'S');
            "Alternative Yeast Nuclear"
        }
        13 => {
            // Ascidian Mitochondrial
            mapping.insert("AGA".into(), 'G');
            mapping.insert("AGG".into(), 'G');
            mapping.insert("ATA".into(), 'M');
            mapping.insert("TGA".into(), 'W');
            "Ascidian Mitochondrial"
        }
        14 => {
            // Alternative Flatworm Mitochondrial
            mapping.insert("AAA".into(), 'N');
            mapping.insert("AGA".into(), 'S');
            mapping.insert("AGG".into(), 'S');
            mapping.insert("TAA".into(), 'Y');
            mapping.insert("TGA".into(), 'W');
            "Alternative Flatworm Mitochondrial"
        }
        15 => {
            // Blepharisma Nuclear
            mapping.insert("TAG".into(), 'Q');
            mapping.insert("TAA".into(), 'Q'); // absent
            "Blepharisma Nuclear"
        }
        16 => {
            // Chlorophycean Mitochondrial
            mapping.insert("TAG".into(), 'L');
            "Chlorophycean Mitochondrial"
        }
        21 => {
            // Trematode Mitochondrial
            mapping.insert("TGA".into(), 'W');
            mapping.insert("ATA".into(), 'M');
            mapping.insert("AGA".into(), 'S');
            mapping.insert("AGG".into(), 'S');
            mapping.insert("AAA".into(), 'N');
            "Trematode Mitochondrial"
        }
        22 => {
            // Scenedesmus obliquus Mitochondrial
            mapping.insert("TCA".into(), '*');
            mapping.insert("TAG".into(), 'L');
            "Scenedesmus obliquus Mitochondrial"
        }
        23 => {
            // Thraustochytrium Mitochondrial
            mapping.insert("TTA".into(), '*');
            mapping.insert("TTG".into(), '*'); // absent
            "Thraustochytrium Mitochondrial"
        }
        24 => {
            // Pterobranchia Mitochondrial
            mapping.insert("AGA".into(), 'S');
            mapping.insert("AGG".into(), 'K');
            mapping.insert("TGA".into(), 'W');
            "Pterobranchia Mitochondrial"
        }
        25 => {
            // Candidate Division SR1/Gracilibacteria
            mapping.insert("TGA".into(), 'G');
            "Candidate Division SR1/Gracilibacteria"
        }
        26 => {
            // Pachysolen tannophilus Nuclear
            mapping.insert("CTG".into(), 'A');
            "Pachysolen tannophilus Nuclear"
        }
        27 => {
            // Karyorelict Nuclear
            mapping.insert("TAA".into(), 'Q');
            mapping.insert("TAG".into(), 'Q');
            mapping.insert("TGA".into(), 'W'); // or 'W'
            "Karyorelict Nuclear"
        }
        28 => {
            // Condylostoma Nuclear
            mapping.insert("TAA".into(), 'Q');
            mapping.insert("TAG".into(), 'Q');
            mapping.insert("TGA".into(), 'W');
            "Condylostoma Nuclear"
        }
        29 => {
            // Mesodinium Nuclear
            mapping.insert("TAA".into(), 'Y');
            mapping.insert("TAG".into(), 'Y');
            "Mesodinium Nuclear"
        }
        30 => {
            // Peritrich Nuclear
            mapping.insert("TAA".into(), 'E');
            mapping.insert("TAG".into(), 'E');
            "Peritrich Nuclear"
        }
        31 => {
            // Blastocrithidia Nuclear
            mapping.insert("TAA".into(), 'E');
            mapping.insert("TAG".into(), 'E');
            mapping.insert("TGA".into(), 'W');
            "Blastocrithidia Nuclear"
        }
        33 => {
            // Cephalodiscidae Mitochondrial
            mapping.insert("AGA".into(), 'S');
            mapping.insert("AGG".into(), 'K');
            mapping.insert("TAA".into(), 'Y');
            mapping.insert("TGA".into(), 'W');
            "Cephalodiscidae Mitochondrial"
        }
        _ => return Err(format!("Unknown NCBI genetic code ID: {gcid}. Available: 1-6, 9-16, 21-31, 33")),
    };

    Ok((name.to_string(), mapping))
}

// ---------------------------------------------------------------
// Standard Genetic Code (NCBI #1)
// ---------------------------------------------------------------

use std::sync::LazyLock;

static STANDARD_CODE: LazyLock<HashMap<String, char>> = LazyLock::new(|| {
    let mut m = HashMap::new();

    // T-start
    m.insert("TTT".into(), 'F'); m.insert("TTC".into(), 'F'); // Phe
    m.insert("TTA".into(), 'L'); m.insert("TTG".into(), 'L'); // Leu
    // C-start
    m.insert("CTT".into(), 'L'); m.insert("CTC".into(), 'L');
    m.insert("CTA".into(), 'L'); m.insert("CTG".into(), 'L'); // Leu
    // A-start
    m.insert("ATT".into(), 'I'); m.insert("ATC".into(), 'I');
    m.insert("ATA".into(), 'I'); // Ile
    m.insert("ATG".into(), 'M'); // Met (start)
    // G-start
    m.insert("GTT".into(), 'V'); m.insert("GTC".into(), 'V');
    m.insert("GTA".into(), 'V'); m.insert("GTG".into(), 'V'); // Val

    // TC
    m.insert("TCT".into(), 'S'); m.insert("TCC".into(), 'S');
    m.insert("TCA".into(), 'S'); m.insert("TCG".into(), 'S'); // Ser
    // CC
    m.insert("CCT".into(), 'P'); m.insert("CCC".into(), 'P');
    m.insert("CCA".into(), 'P'); m.insert("CCG".into(), 'P'); // Pro
    // AC
    m.insert("ACT".into(), 'T'); m.insert("ACC".into(), 'T');
    m.insert("ACA".into(), 'T'); m.insert("ACG".into(), 'T'); // Thr
    // GC
    m.insert("GCT".into(), 'A'); m.insert("GCC".into(), 'A');
    m.insert("GCA".into(), 'A'); m.insert("GCG".into(), 'A'); // Ala

    // TA
    m.insert("TAT".into(), 'Y'); m.insert("TAC".into(), 'Y'); // Tyr
    m.insert("TAA".into(), '*'); m.insert("TAG".into(), '*'); // Stop
    // CA
    m.insert("CAT".into(), 'H'); m.insert("CAC".into(), 'H'); // His
    m.insert("CAA".into(), 'Q'); m.insert("CAG".into(), 'Q'); // Gln
    // AA
    m.insert("AAT".into(), 'N'); m.insert("AAC".into(), 'N'); // Asn
    m.insert("AAA".into(), 'K'); m.insert("AAG".into(), 'K'); // Lys
    // GA
    m.insert("GAT".into(), 'D'); m.insert("GAC".into(), 'D'); // Asp
    m.insert("GAA".into(), 'E'); m.insert("GAG".into(), 'E'); // Glu

    // TG
    m.insert("TGT".into(), 'C'); m.insert("TGC".into(), 'C'); // Cys
    m.insert("TGA".into(), '*'); // Stop (Opal)
    m.insert("TGG".into(), 'W'); // Trp
    // CG
    m.insert("CGT".into(), 'R'); m.insert("CGC".into(), 'R');
    m.insert("CGA".into(), 'R'); m.insert("CGG".into(), 'R'); // Arg
    // AG
    m.insert("AGT".into(), 'S'); m.insert("AGC".into(), 'S'); // Ser
    m.insert("AGA".into(), 'R'); m.insert("AGG".into(), 'R'); // Arg
    // GG
    m.insert("GGT".into(), 'G'); m.insert("GGC".into(), 'G');
    m.insert("GGA".into(), 'G'); m.insert("GGG".into(), 'G'); // Gly

    m
});

// ================================================================
// Tests
// ================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_code() {
        let ct = CodonTable::standard();
        assert_eq!(ct.gcid, "1");
        assert_eq!(ct.get_aa("ATG").unwrap(), 'M');
        assert_eq!(ct.get_aa("TGG").unwrap(), 'W');
        assert_eq!(ct.get_aa("TAA").unwrap(), '*');
        assert!(ct.is_start("ATG"));
        assert!(ct.is_stop("TAA"));
        assert!(ct.is_stop("TAG"));
        assert!(ct.is_stop("TGA"));

        // Leucine has 6 codons
        let leu = ct.synonymous_codons('L');
        assert_eq!(leu.len(), 6);
        assert!(leu.contains(&"TTA".to_string()));
        assert!(leu.contains(&"TTG".to_string()));
        assert!(leu.contains(&"CTT".to_string()));

        // Serine has 6 codons
        let ser = ct.synonymous_codons('S');
        assert_eq!(ser.len(), 6);

        // Methionine has 1 codon
        let met = ct.synonymous_codons('M');
        assert_eq!(met.len(), 1);
    }

    #[test]
    fn test_subfamilies() {
        let ct = CodonTable::standard();
        let groups = ct.subfamily_groups();

        // Leucine should have 2 subfamilies
        assert!(groups.contains_key("Leu_TT"));
        assert!(groups.contains_key("Leu_CT"));
        assert_eq!(groups["Leu_TT"].len(), 2); // TTA, TTG
        assert_eq!(groups["Leu_CT"].len(), 4); // CTT, CTC, CTA, CTG

        // Serine should have 2 subfamilies
        assert!(groups.contains_key("Ser_TC"));
        assert!(groups.contains_key("Ser_AG"));
    }

    #[test]
    fn test_vertebrate_mitochondrial() {
        let ct = CodonTable::from_ncbi_id("2").unwrap();
        // AGA/AGG are stop in vertebrate mitochondrial
        assert_eq!(ct.get_aa("AGA").unwrap(), '*');
        assert_eq!(ct.get_aa("AGG").unwrap(), '*');
        // ATA is Met
        assert_eq!(ct.get_aa("ATA").unwrap(), 'M');
        // TGA is Trp
        assert_eq!(ct.get_aa("TGA").unwrap(), 'W');
    }

    #[test]
    fn test_list_codes() {
        let codes = CodonTable::list_available();
        assert!(codes.iter().any(|(id, _)| id == "1"));
        assert!(codes.iter().any(|(id, _)| id == "2"));
    }
}
