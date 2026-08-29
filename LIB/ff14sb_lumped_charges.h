#ifndef FF14SB_LUMPED_CHARGES_H
#define FF14SB_LUMPED_CHARGES_H
// GENERATED FILE - DO NOT EDIT BY HAND.
// Source: openmm/app/data/amber14/protein.ff14SB.xml (archived as an artifact alongside the
//   table it produced). Every value below is machine-extracted; none was typed by hand.
// Transform: HYDROGEN LUMPING - each hydrogen's partial charge is folded into the heavy atom
//   it is bonded to, producing a united-atom set that matches this engine's heavy-atom-only
//   representation and the ligand-side Gasteiger preparation (both sides then use ONE
//   convention, which is the defect this file exists to fix).
// Gates, both verified at generation time:
//   integrality  - every ff14SB template sums to an integer charge
//   conservation - lumped sum == all-atom sum for all 70 templates, max residual < 1e-6
// Coverage: 24 standard residues + 23 N-terminal + 23 C-terminal templates.
// HIS: the PDB writes "HIS" while ff14SB requires a tautomer/protonation choice. "HIS" is
//   aliased here to HIE (neutral, NE2-protonated). HID, HIE and HIP are ALSO emitted so an
//   explicitly-named residue gets its own values. This alias is a DECLARED MODELLING CHOICE,
//   not a measurement - ND1/NE2 differ by up to 0.79 e between the three states.
// Metals: NOT covered here. METAL_ION_CHARGES in assign_formal_charges.h keeps its formal
//   integer charges, which is itself a convention mismatch against these partials.
namespace ff14sb_lumped {
struct Entry { const char* res_name; const char* atom_name; float charge; };
static constexpr Entry FF14SB_LUMPED_CHARGES[] = {
    { "ALA", " C  ", +0.5973f },   // STD
    { "ALA", " CA ", +0.1160f },   // STD
    { "ALA", " CB ", -0.0016f },   // STD
    { "ALA", " N  ", -0.1438f },   // STD
    { "ALA", " O  ", -0.5679f },   // STD
    { "ARG", " C  ", +0.7341f },   // STD
    { "ARG", " CA ", -0.1077f },   // STD
    { "ARG", " CB ", +0.0647f },   // STD
    { "ARG", " CD ", +0.1860f },   // STD
    { "ARG", " CG ", +0.0960f },   // STD
    { "ARG", " CZ ", +0.8076f },   // STD
    { "ARG", " N  ", -0.0732f },   // STD
    { "ARG", " NE ", -0.1839f },   // STD
    { "ARG", " NH1", +0.0329f },   // STD
    { "ARG", " NH2", +0.0329f },   // STD
    { "ARG", " O  ", -0.5894f },   // STD
    { "ASN", " C  ", +0.5973f },   // STD
    { "ASN", " CA ", +0.1191f },   // STD
    { "ASN", " CB ", -0.0447f },   // STD
    { "ASN", " CG ", +0.7130f },   // STD
    { "ASN", " N  ", -0.1438f },   // STD
    { "ASN", " ND2", -0.0799f },   // STD
    { "ASN", " O  ", -0.5679f },   // STD
    { "ASN", " OD1", -0.5931f },   // STD
    { "ASP", " C  ", +0.5366f },   // STD
    { "ASP", " CA ", +0.1261f },   // STD
    { "ASP", " CB ", -0.0547f },   // STD
    { "ASP", " CG ", +0.7994f },   // STD
    { "ASP", " N  ", -0.2227f },   // STD
    { "ASP", " O  ", -0.5819f },   // STD
    { "ASP", " OD1", -0.8014f },   // STD
    { "ASP", " OD2", -0.8014f },   // STD
    { "CYS", " C  ", +0.5973f },   // STD
    { "CYS", " CA ", +0.1337f },   // STD
    { "CYS", " CB ", +0.0993f },   // STD
    { "CYS", " N  ", -0.1438f },   // STD
    { "CYS", " O  ", -0.5679f },   // STD
    { "CYS", " SG ", -0.1186f },   // STD
    { "CYX", " C  ", +0.5973f },   // STD
    { "CYX", " CA ", +0.1195f },   // STD
    { "CYX", " CB ", +0.1030f },   // STD
    { "CYX", " N  ", -0.1438f },   // STD
    { "CYX", " O  ", -0.5679f },   // STD
    { "CYX", " SG ", -0.1081f },   // STD
    { "CYM", " C  ", +0.5973f },   // STD
    { "CYM", " CA ", +0.0157f },   // STD
    { "CYM", " CB ", -0.0169f },   // STD
    { "CYM", " N  ", -0.1438f },   // STD
    { "CYM", " O  ", -0.5679f },   // STD
    { "CYM", " SG ", -0.8844f },   // STD
    { "GLN", " C  ", +0.5973f },   // STD
    { "GLN", " CA ", +0.0819f },   // STD
    { "GLN", " CB ", +0.0306f },   // STD
    { "GLN", " CD ", +0.6951f },   // STD
    { "GLN", " CG ", +0.0059f },   // STD
    { "GLN", " N  ", -0.1438f },   // STD
    { "GLN", " NE2", -0.0905f },   // STD
    { "GLN", " O  ", -0.5679f },   // STD
    { "GLN", " OE1", -0.6086f },   // STD
    { "GLU", " C  ", +0.5366f },   // STD
    { "GLU", " CA ", +0.1502f },   // STD
    { "GLU", " CB ", +0.0214f },   // STD
    { "GLU", " CD ", +0.8054f },   // STD
    { "GLU", " CG ", -0.0714f },   // STD
    { "GLU", " N  ", -0.2227f },   // STD
    { "GLU", " O  ", -0.5819f },   // STD
    { "GLU", " OE1", -0.8188f },   // STD
    { "GLU", " OE2", -0.8188f },   // STD
    { "GLY", " C  ", +0.5973f },   // STD
    { "GLY", " CA ", +0.1144f },   // STD
    { "GLY", " N  ", -0.1438f },   // STD
    { "GLY", " O  ", -0.5679f },   // STD
    { "HID", " C  ", +0.5973f },   // STD
    { "HID", " CA ", +0.1069f },   // STD
    { "HID", " CB ", +0.0342f },   // STD
    { "HID", " CD2", +0.2439f },   // STD
    { "HID", " CE1", +0.3449f },   // STD
    { "HID", " CG ", -0.0266f },   // STD
    { "HID", " N  ", -0.1438f },   // STD
    { "HID", " ND1", -0.0162f },   // STD
    { "HID", " NE2", -0.5727f },   // STD
    { "HID", " O  ", -0.5679f },   // STD
    { "HIE", " C  ", +0.5973f },   // STD
    { "HIE", " CA ", +0.0779f },   // STD
    { "HIE", " CB ", +0.0660f },   // STD
    { "HIE", " CD2", -0.0345f },   // STD
    { "HIE", " CE1", +0.3070f },   // STD
    { "HIE", " CG ", +0.1868f },   // STD
    { "HIE", " N  ", -0.1438f },   // STD
    { "HIE", " ND1", -0.5432f },   // STD
    { "HIE", " NE2", +0.0544f },   // STD
    { "HIE", " O  ", -0.5679f },   // STD
    { "HIP", " C  ", +0.7341f },   // STD
    { "HIP", " CA ", -0.0142f },   // STD
    { "HIP", " CB ", +0.1206f },   // STD
    { "HIP", " CD2", +0.1176f },   // STD
    { "HIP", " CE1", +0.2511f },   // STD
    { "HIP", " CG ", -0.0012f },   // STD
    { "HIP", " N  ", -0.0732f },   // STD
    { "HIP", " ND1", +0.2353f },   // STD
    { "HIP", " NE2", +0.2193f },   // STD
    { "HIP", " O  ", -0.5894f },   // STD
    { "ILE", " C  ", +0.5973f },   // STD
    { "ILE", " CA ", +0.0272f },   // STD
    { "ILE", " CB ", +0.1490f },   // STD
    { "ILE", " CD1", -0.0102f },   // STD
    { "ILE", " CG1", +0.0042f },   // STD
    { "ILE", " CG2", -0.0558f },   // STD
    { "ILE", " N  ", -0.1438f },   // STD
    { "ILE", " O  ", -0.5679f },   // STD
    { "LEU", " C  ", +0.5973f },   // STD
    { "LEU", " CA ", +0.0404f },   // STD
    { "LEU", " CB ", -0.0188f },   // STD
    { "LEU", " CD1", -0.1121f },   // STD
    { "LEU", " CD2", -0.1121f },   // STD
    { "LEU", " CG ", +0.3170f },   // STD
    { "LEU", " N  ", -0.1438f },   // STD
    { "LEU", " O  ", -0.5679f },   // STD
    { "LYS", " C  ", +0.7341f },   // STD
    { "LYS", " CA ", -0.0974f },   // STD
    { "LYS", " CB ", +0.0630f },   // STD
    { "LYS", " CD ", +0.0763f },   // STD
    { "LYS", " CE ", +0.2127f },   // STD
    { "LYS", " CG ", +0.0393f },   // STD
    { "LYS", " N  ", -0.0732f },   // STD
    { "LYS", " NZ ", +0.6346f },   // STD
    { "LYS", " O  ", -0.5894f },   // STD
    { "MET", " C  ", +0.5973f },   // STD
    { "MET", " CA ", +0.0643f },   // STD
    { "MET", " CB ", +0.0824f },   // STD
    { "MET", " CE ", +0.1516f },   // STD
    { "MET", " CG ", +0.0898f },   // STD
    { "MET", " N  ", -0.1438f },   // STD
    { "MET", " O  ", -0.5679f },   // STD
    { "MET", " SD ", -0.2737f },   // STD
    { "PHE", " C  ", +0.5973f },   // STD
    { "PHE", " CA ", +0.0954f },   // STD
    { "PHE", " CB ", +0.0247f },   // STD
    { "PHE", " CD1", +0.0074f },   // STD
    { "PHE", " CD2", +0.0074f },   // STD
    { "PHE", " CE1", -0.0274f },   // STD
    { "PHE", " CE2", -0.0274f },   // STD
    { "PHE", " CG ", +0.0118f },   // STD
    { "PHE", " CZ ", +0.0225f },   // STD
    { "PHE", " N  ", -0.1438f },   // STD
    { "PHE", " O  ", -0.5679f },   // STD
    { "PRO", " C  ", +0.5896f },   // STD
    { "PRO", " CA ", +0.0375f },   // STD
    { "PRO", " CB ", +0.0436f },   // STD
    { "PRO", " CD ", +0.0974f },   // STD
    { "PRO", " CG ", +0.0615f },   // STD
    { "PRO", " N  ", -0.2548f },   // STD
    { "PRO", " O  ", -0.5748f },   // STD
    { "SER", " C  ", +0.5973f },   // STD
    { "SER", " CA ", +0.0594f },   // STD
    { "SER", " CB ", +0.2821f },   // STD
    { "SER", " N  ", -0.1438f },   // STD
    { "SER", " O  ", -0.5679f },   // STD
    { "SER", " OG ", -0.2271f },   // STD
    { "THR", " C  ", +0.5973f },   // STD
    { "THR", " CA ", +0.0618f },   // STD
    { "THR", " CB ", +0.3697f },   // STD
    { "THR", " CG2", -0.0512f },   // STD
    { "THR", " N  ", -0.1438f },   // STD
    { "THR", " O  ", -0.5679f },   // STD
    { "THR", " OG1", -0.2659f },   // STD
    { "TRP", " C  ", +0.5973f },   // STD
    { "TRP", " CA ", +0.0848f },   // STD
    { "TRP", " CB ", +0.0628f },   // STD
    { "TRP", " CD1", +0.0424f },   // STD
    { "TRP", " CD2", +0.1243f },   // STD
    { "TRP", " CE2", +0.1380f },   // STD
    { "TRP", " CE3", -0.0687f },   // STD
    { "TRP", " CG ", -0.1415f },   // STD
    { "TRP", " CH2", +0.0283f },   // STD
    { "TRP", " CZ2", -0.1029f },   // STD
    { "TRP", " CZ3", -0.0525f },   // STD
    { "TRP", " N  ", -0.1438f },   // STD
    { "TRP", " NE1", -0.0006f },   // STD
    { "TRP", " O  ", -0.5679f },   // STD
    { "TYR", " C  ", +0.5973f },   // STD
    { "TYR", " CA ", +0.0862f },   // STD
    { "TYR", " CB ", +0.0438f },   // STD
    { "TYR", " CD1", -0.0207f },   // STD
    { "TYR", " CD2", -0.0207f },   // STD
    { "TYR", " CE1", -0.0685f },   // STD
    { "TYR", " CE2", -0.0685f },   // STD
    { "TYR", " CG ", -0.0011f },   // STD
    { "TYR", " CZ ", +0.3226f },   // STD
    { "TYR", " N  ", -0.1438f },   // STD
    { "TYR", " O  ", -0.5679f },   // STD
    { "TYR", " OH ", -0.1587f },   // STD
    { "VAL", " C  ", +0.5973f },   // STD
    { "VAL", " CA ", +0.0094f },   // STD
    { "VAL", " CB ", +0.2688f },   // STD
    { "VAL", " CG1", -0.0819f },   // STD
    { "VAL", " CG2", -0.0819f },   // STD
    { "VAL", " N  ", -0.1438f },   // STD
    { "VAL", " O  ", -0.5679f },   // STD
    { "NALA", " C  ", +0.6163f },   // NTERM
    { "NALA", " CA ", +0.1851f },   // NTERM
    { "NALA", " CB ", +0.0303f },   // NTERM
    { "NALA", " N  ", +0.7405f },   // NTERM
    { "NALA", " O  ", -0.5722f },   // NTERM
    { "NARG", " C  ", +0.7214f },   // NTERM
    { "NARG", " CA ", +0.1019f },   // NTERM
    { "NARG", " CB ", +0.0570f },   // NTERM
    { "NARG", " CD ", +0.1989f },   // NTERM
    { "NARG", " CG ", +0.0854f },   // NTERM
    { "NARG", " CZ ", +0.8281f },   // NTERM
    { "NARG", " N  ", +0.7554f },   // NTERM
    { "NARG", " NE ", -0.2058f },   // NTERM
    { "NARG", " NH1", +0.0295f },   // NTERM
    { "NARG", " NH2", +0.0295f },   // NTERM
    { "NARG", " O  ", -0.6013f },   // NTERM
    { "NASN", " C  ", +0.6163f },   // NTERM
    { "NASN", " CA ", +0.1599f },   // NTERM
    { "NASN", " CB ", +0.0747f },   // NTERM
    { "NASN", " CG ", +0.5833f },   // NTERM
    { "NASN", " N  ", +0.7564f },   // NTERM
    { "NASN", " ND2", -0.0440f },   // NTERM
    { "NASN", " O  ", -0.5722f },   // NTERM
    { "NASN", " OD1", -0.5744f },   // NTERM
    { "NASP", " C  ", +0.5621f },   // NTERM
    { "NASP", " CA ", +0.1433f },   // NTERM
    { "NASP", " CB ", -0.0573f },   // NTERM
    { "NASP", " CG ", +0.8194f },   // NTERM
    { "NASP", " N  ", +0.7382f },   // NTERM
    { "NASP", " O  ", -0.5889f },   // NTERM
    { "NASP", " OD1", -0.8084f },   // NTERM
    { "NASP", " OD2", -0.8084f },   // NTERM
    { "NCYS", " C  ", +0.6123f },   // NTERM
    { "NCYS", " CA ", +0.2338f },   // NTERM
    { "NCYS", " CB ", +0.1181f },   // NTERM
    { "NCYS", " N  ", +0.7394f },   // NTERM
    { "NCYS", " O  ", -0.5713f },   // NTERM
    { "NCYS", " SG ", -0.1323f },   // NTERM
    { "NCYX", " C  ", +0.6123f },   // NTERM
    { "NCYX", " CA ", +0.1977f },   // NTERM
    { "NCYX", " CB ", +0.1083f },   // NTERM
    { "NCYX", " N  ", +0.7514f },   // NTERM
    { "NCYX", " O  ", -0.5713f },   // NTERM
    { "NCYX", " SG ", -0.0984f },   // NTERM
    { "NGLN", " C  ", +0.6123f },   // NTERM
    { "NGLN", " CA ", +0.1551f },   // NTERM
    { "NGLN", " CB ", +0.0751f },   // NTERM
    { "NGLN", " CD ", +0.7354f },   // NTERM
    { "NGLN", " CG ", -0.0241f },   // NTERM
    { "NGLN", " N  ", +0.7481f },   // NTERM
    { "NGLN", " NE2", -0.1173f },   // NTERM
    { "NGLN", " O  ", -0.5713f },   // NTERM
    { "NGLN", " OE1", -0.6133f },   // NTERM
    { "NGLU", " C  ", +0.5621f },   // NTERM
    { "NGLU", " CA ", +0.1790f },   // NTERM
    { "NGLU", " CB ", +0.0445f },   // NTERM
    { "NGLU", " CD ", +0.8087f },   // NTERM
    { "NGLU", " CG ", -0.0866f },   // NTERM
    { "NGLU", " N  ", +0.7190f },   // NTERM
    { "NGLU", " O  ", -0.5889f },   // NTERM
    { "NGLU", " OE1", -0.8189f },   // NTERM
    { "NGLU", " OE2", -0.8189f },   // NTERM
    { "NGLY", " C  ", +0.6163f },   // NTERM
    { "NGLY", " CA ", +0.1690f },   // NTERM
    { "NGLY", " N  ", +0.7869f },   // NTERM
    { "NGLY", " O  ", -0.5722f },   // NTERM
    { "NHID", " C  ", +0.6123f },   // NTERM
    { "NHID", " CA ", +0.1922f },   // NTERM
    { "NHID", " CB ", +0.0677f },   // NTERM
    { "NHID", " CD2", +0.2345f },   // NTERM
    { "NHID", " CE1", +0.3512f },   // NTERM
    { "NHID", " CG ", -0.0399f },   // NTERM
    { "NHID", " N  ", +0.7431f },   // NTERM
    { "NHID", " ND1", -0.0187f },   // NTERM
    { "NHID", " NE2", -0.5711f },   // NTERM
    { "NHID", " O  ", -0.5713f },   // NTERM
    { "NHIE", " C  ", +0.6123f },   // NTERM
    { "NHIE", " CA ", +0.1616f },   // NTERM
    { "NHIE", " CB ", +0.0935f },   // NTERM
    { "NHIE", " CD2", -0.0386f },   // NTERM
    { "NHIE", " CE1", +0.3201f },   // NTERM
    { "NHIE", " CG ", +0.1740f },   // NTERM
    { "NHIE", " N  ", +0.7520f },   // NTERM
    { "NHIE", " ND1", -0.5579f },   // NTERM
    { "NHIE", " NE2", +0.0543f },   // NTERM
    { "NHIE", " O  ", -0.5713f },   // NTERM
    { "NHIP", " C  ", +0.7214f },   // NTERM
    { "NHIP", " CA ", +0.1628f },   // NTERM
    { "NHIP", " CB ", +0.1546f },   // NTERM
    { "NHIP", " CD2", +0.1062f },   // NTERM
    { "NHIP", " CE1", +0.2634f },   // NTERM
    { "NHIP", " CG ", -0.0236f },   // NTERM
    { "NHIP", " N  ", +0.7672f },   // NTERM
    { "NHIP", " ND1", +0.2311f },   // NTERM
    { "NHIP", " NE2", +0.2182f },   // NTERM
    { "NHIP", " O  ", -0.6013f },   // NTERM
    { "NILE", " C  ", +0.6123f },   // NTERM
    { "NILE", " CA ", +0.1288f },   // NTERM
    { "NILE", " CB ", +0.2098f },   // NTERM
    { "NILE", " CD1", -0.0230f },   // NTERM
    { "NILE", " CG1", +0.0015f },   // NTERM
    { "NILE", " CG2", -0.0879f },   // NTERM
    { "NILE", " N  ", +0.7298f },   // NTERM
    { "NILE", " O  ", -0.5713f },   // NTERM
    { "NLEU", " C  ", +0.6123f },   // NTERM
    { "NLEU", " CA ", +0.1157f },   // NTERM
    { "NLEU", " CB ", +0.0268f },   // NTERM
    { "NLEU", " CD1", -0.1166f },   // NTERM
    { "NLEU", " CD2", -0.1164f },   // NTERM
    { "NLEU", " CG ", +0.3041f },   // NTERM
    { "NLEU", " N  ", +0.7454f },   // NTERM
    { "NLEU", " O  ", -0.5713f },   // NTERM
    { "NLYS", " C  ", +0.7214f },   // NTERM
    { "NLYS", " CA ", +0.1165f },   // NTERM
    { "NLYS", " CB ", +0.0778f },   // NTERM
    { "NLYS", " CD ", +0.0658f },   // NTERM
    { "NLYS", " CE ", +0.2161f },   // NTERM
    { "NLYS", " CG ", +0.0194f },   // NTERM
    { "NLYS", " N  ", +0.7461f },   // NTERM
    { "NLYS", " NZ ", +0.6382f },   // NTERM
    { "NLYS", " O  ", -0.6013f },   // NTERM
    { "NMET", " C  ", +0.6123f },   // NTERM
    { "NMET", " CA ", +0.1337f },   // NTERM
    { "NMET", " CB ", +0.1115f },   // NTERM
    { "NMET", " CE ", +0.1450f },   // NTERM
    { "NMET", " CG ", +0.0918f },   // NTERM
    { "NMET", " N  ", +0.7544f },   // NTERM
    { "NMET", " O  ", -0.5713f },   // NTERM
    { "NMET", " SD ", -0.2774f },   // NTERM
    { "NPHE", " C  ", +0.6123f },   // NTERM
    { "NPHE", " CA ", +0.1774f },   // NTERM
    { "NPHE", " CB ", +0.0538f },   // NTERM
    { "NPHE", " CD1", -0.0018f },   // NTERM
    { "NPHE", " CD2", -0.0017f },   // NTERM
    { "NPHE", " CE1", -0.0169f },   // NTERM
    { "NPHE", " CE2", -0.0170f },   // NTERM
    { "NPHE", " CG ", +0.0031f },   // NTERM
    { "NPHE", " CZ ", +0.0121f },   // NTERM
    { "NPHE", " N  ", +0.7500f },   // NTERM
    { "NPHE", " O  ", -0.5713f },   // NTERM
    { "NPRO", " C  ", +0.5260f },   // NTERM
    { "NPRO", " CA ", +0.2000f },   // NTERM
    { "NPRO", " CB ", +0.0850f },   // NTERM
    { "NPRO", " CD ", +0.1880f },   // NTERM
    { "NPRO", " CG ", +0.0790f },   // NTERM
    { "NPRO", " N  ", +0.4220f },   // NTERM
    { "NPRO", " O  ", -0.5000f },   // NTERM
    { "NSER", " C  ", +0.6163f },   // NTERM
    { "NSER", " CA ", +0.1349f },   // NTERM
    { "NSER", " CB ", +0.3142f },   // NTERM
    { "NSER", " N  ", +0.7543f },   // NTERM
    { "NSER", " O  ", -0.5722f },   // NTERM
    { "NSER", " OG ", -0.2475f },   // NTERM
    { "NTHR", " C  ", +0.6163f },   // NTERM
    { "NTHR", " CA ", +0.1121f },   // NTERM
    { "NTHR", " CB ", +0.4191f },   // NTERM
    { "NTHR", " CG2", -0.0673f },   // NTERM
    { "NTHR", " N  ", +0.7614f },   // NTERM
    { "NTHR", " O  ", -0.5722f },   // NTERM
    { "NTHR", " OG1", -0.2694f },   // NTERM
    { "NTRP", " C  ", +0.6123f },   // NTERM
    { "NTRP", " CA ", +0.1583f },   // NTERM
    { "NTRP", " CB ", +0.0987f },   // NTERM
    { "NTRP", " CD1", +0.0407f },   // NTERM
    { "NTRP", " CD2", +0.1132f },   // NTERM
    { "NTRP", " CE2", +0.1575f },   // NTERM
    { "NTRP", " CE3", -0.0619f },   // NTERM
    { "NTRP", " CG ", -0.1654f },   // NTERM
    { "NTRP", " CH2", +0.0331f },   // NTERM
    { "NTRP", " CZ2", -0.1121f },   // NTERM
    { "NTRP", " CZ3", -0.0576f },   // NTERM
    { "NTRP", " N  ", +0.7577f },   // NTERM
    { "NTRP", " NE1", -0.0032f },   // NTERM
    { "NTRP", " O  ", -0.5713f },   // NTERM
    { "NTYR", " C  ", +0.6123f },   // NTERM
    { "NTYR", " CA ", +0.1553f },   // NTERM
    { "NTYR", " CB ", +0.0863f },   // NTERM
    { "NTYR", " CD1", -0.0282f },   // NTERM
    { "NTYR", " CD2", -0.0282f },   // NTERM
    { "NTYR", " CE1", -0.0589f },   // NTERM
    { "NTYR", " CE2", -0.0589f },   // NTERM
    { "NTYR", " CG ", -0.0205f },   // NTERM
    { "NTYR", " CZ ", +0.3139f },   // NTERM
    { "NTYR", " N  ", +0.7559f },   // NTERM
    { "NTYR", " O  ", -0.5713f },   // NTERM
    { "NTYR", " OH ", -0.1577f },   // NTERM
    { "NVAL", " C  ", +0.6163f },   // NTERM
    { "NVAL", " CA ", +0.1039f },   // NTERM
    { "NVAL", " CB ", +0.2975f },   // NTERM
    { "NVAL", " CG1", -0.0924f },   // NTERM
    { "NVAL", " CG2", -0.0924f },   // NTERM
    { "NVAL", " N  ", +0.7393f },   // NTERM
    { "NVAL", " O  ", -0.5722f },   // NTERM
    { "CALA", " C  ", +0.7731f },   // CTERM
    { "CALA", " CA ", -0.0680f },   // CTERM
    { "CALA", " CB ", +0.0199f },   // CTERM
    { "CALA", " N  ", -0.1140f },   // CTERM
    { "CALA", " O  ", -0.8055f },   // CTERM
    { "CALA", " OXT", -0.8055f },   // CTERM
    { "CARG", " C  ", +0.8557f },   // CTERM
    { "CARG", " CA ", -0.1621f },   // CTERM
    { "CARG", " CB ", +0.0368f },   // CTERM
    { "CARG", " CD ", +0.2050f },   // CTERM
    { "CARG", " CG ", +0.1114f },   // CTERM
    { "CARG", " CZ ", +0.8368f },   // CTERM
    { "CARG", " N  ", -0.0717f },   // CTERM
    { "CARG", " NE ", -0.2085f },   // CTERM
    { "CARG", " NH1", +0.0249f },   // CTERM
    { "CARG", " NH2", +0.0249f },   // CTERM
    { "CARG", " O  ", -0.8266f },   // CTERM
    { "CARG", " OXT", -0.8266f },   // CTERM
    { "CASN", " C  ", +0.8050f },   // CTERM
    { "CASN", " CA ", -0.0722f },   // CTERM
    { "CASN", " CB ", -0.0253f },   // CTERM
    { "CASN", " CG ", +0.7153f },   // CTERM
    { "CASN", " N  ", -0.1140f },   // CTERM
    { "CASN", " ND2", -0.0784f },   // CTERM
    { "CASN", " O  ", -0.8147f },   // CTERM
    { "CASN", " OD1", -0.6010f },   // CTERM
    { "CASN", " OXT", -0.8147f },   // CTERM
    { "CASP", " C  ", +0.7256f },   // CTERM
    { "CASP", " CA ", -0.0771f },   // CTERM
    { "CASP", " CB ", -0.1101f },   // CTERM
    { "CASP", " CG ", +0.8851f },   // CTERM
    { "CASP", " N  ", -0.2137f },   // CTERM
    { "CASP", " O  ", -0.7887f },   // CTERM
    { "CASP", " OD1", -0.8162f },   // CTERM
    { "CASP", " OD2", -0.8162f },   // CTERM
    { "CASP", " OXT", -0.7887f },   // CTERM
    { "CCYS", " C  ", +0.7497f },   // CTERM
    { "CCYS", " CA ", -0.0239f },   // CTERM
    { "CCYS", " CB ", +0.0878f },   // CTERM
    { "CCYS", " N  ", -0.1140f },   // CTERM
    { "CCYS", " O  ", -0.7981f },   // CTERM
    { "CCYS", " OXT", -0.7981f },   // CTERM
    { "CCYS", " SG ", -0.1034f },   // CTERM
    { "CCYX", " C  ", +0.7618f },   // CTERM
    { "CCYX", " CA ", -0.0380f },   // CTERM
    { "CCYX", " CB ", +0.0513f },   // CTERM
    { "CCYX", " N  ", -0.1140f },   // CTERM
    { "CCYX", " O  ", -0.8041f },   // CTERM
    { "CCYX", " OXT", -0.8041f },   // CTERM
    { "CCYX", " SG ", -0.0529f },   // CTERM
    { "CGLN", " C  ", +0.7775f },   // CTERM
    { "CGLN", " CA ", -0.1016f },   // CTERM
    { "CGLN", " CB ", +0.0240f },   // CTERM
    { "CGLN", " CD ", +0.7093f },   // CTERM
    { "CGLN", " CG ", +0.0196f },   // CTERM
    { "CGLN", " N  ", -0.1140f },   // CTERM
    { "CGLN", " NE2", -0.0966f },   // CTERM
    { "CGLN", " O  ", -0.8042f },   // CTERM
    { "CGLN", " OE1", -0.6098f },   // CTERM
    { "CGLN", " OXT", -0.8042f },   // CTERM
    { "CGLU", " C  ", +0.7420f },   // CTERM
    { "CGLU", " CA ", -0.0660f },   // CTERM
    { "CGLU", " CB ", -0.0085f },   // CTERM
    { "CGLU", " CD ", +0.8183f },   // CTERM
    { "CGLU", " CG ", -0.0421f },   // CTERM
    { "CGLU", " N  ", -0.2137f },   // CTERM
    { "CGLU", " O  ", -0.7930f },   // CTERM
    { "CGLU", " OE1", -0.8220f },   // CTERM
    { "CGLU", " OE2", -0.8220f },   // CTERM
    { "CGLU", " OXT", -0.7930f },   // CTERM
    { "CGLY", " C  ", +0.7231f },   // CTERM
    { "CGLY", " CA ", -0.0381f },   // CTERM
    { "CGLY", " N  ", -0.1140f },   // CTERM
    { "CGLY", " O  ", -0.7855f },   // CTERM
    { "CGLY", " OXT", -0.7855f },   // CTERM
    { "CHID", " C  ", +0.7615f },   // CTERM
    { "CHID", " CA ", -0.0639f },   // CTERM
    { "CHID", " CB ", +0.0084f },   // CTERM
    { "CHID", " CD2", +0.2242f },   // CTERM
    { "CHID", " CE1", +0.3343f },   // CTERM
    { "CHID", " CG ", +0.0293f },   // CTERM
    { "CHID", " N  ", -0.1140f },   // CTERM
    { "CHID", " ND1", -0.0137f },   // CTERM
    { "CHID", " NE2", -0.5629f },   // CTERM
    { "CHID", " O  ", -0.8016f },   // CTERM
    { "CHID", " OXT", -0.8016f },   // CTERM
    { "CHIE", " C  ", +0.7916f },   // CTERM
    { "CHIE", " CA ", -0.1049f },   // CTERM
    { "CHIE", " CB ", +0.0172f },   // CTERM
    { "CHIE", " CD2", -0.0631f },   // CTERM
    { "CHIE", " CE1", +0.3006f },   // CTERM
    { "CHIE", " CG ", +0.2724f },   // CTERM
    { "CHIE", " N  ", -0.1140f },   // CTERM
    { "CHIE", " ND1", -0.5517f },   // CTERM
    { "CHIE", " NE2", +0.0649f },   // CTERM
    { "CHIE", " O  ", -0.8065f },   // CTERM
    { "CHIE", " OXT", -0.8065f },   // CTERM
    { "CHIP", " C  ", +0.8032f },   // CTERM
    { "CHIP", " CA ", -0.0330f },   // CTERM
    { "CHIP", " CB ", +0.0936f },   // CTERM
    { "CHIP", " CD2", +0.1080f },   // CTERM
    { "CHIP", " CE1", +0.2443f },   // CTERM
    { "CHIP", " CG ", +0.0298f },   // CTERM
    { "CHIP", " N  ", -0.0717f },   // CTERM
    { "CHIP", " ND1", +0.2382f },   // CTERM
    { "CHIP", " NE2", +0.2230f },   // CTERM
    { "CHIP", " O  ", -0.8177f },   // CTERM
    { "CHIP", " OXT", -0.8177f },   // CTERM
    { "CILE", " C  ", +0.8343f },   // CTERM
    { "CILE", " CA ", -0.1725f },   // CTERM
    { "CILE", " CB ", +0.1129f },   // CTERM
    { "CILE", " CD1", -0.0111f },   // CTERM
    { "CILE", " CG1", +0.0319f },   // CTERM
    { "CILE", " CG2", -0.0435f },   // CTERM
    { "CILE", " N  ", -0.1140f },   // CTERM
    { "CILE", " O  ", -0.8190f },   // CTERM
    { "CILE", " OXT", -0.8190f },   // CTERM
    { "CLEU", " C  ", +0.8326f },   // CTERM
    { "CLEU", " CA ", -0.1501f },   // CTERM
    { "CLEU", " CB ", -0.0521f },   // CTERM
    { "CLEU", " CD1", -0.1049f },   // CTERM
    { "CLEU", " CD2", -0.1049f },   // CTERM
    { "CLEU", " CG ", +0.3332f },   // CTERM
    { "CLEU", " N  ", -0.1140f },   // CTERM
    { "CLEU", " O  ", -0.8199f },   // CTERM
    { "CLEU", " OXT", -0.8199f },   // CTERM
    { "CLYS", " C  ", +0.8488f },   // CTERM
    { "CLYS", " CA ", -0.1465f },   // CTERM
    { "CLYS", " CB ", +0.0426f },   // CTERM
    { "CLYS", " CD ", +0.0830f },   // CTERM
    { "CLYS", " CE ", +0.2066f },   // CTERM
    { "CLYS", " CG ", +0.0495f },   // CTERM
    { "CLYS", " N  ", -0.0717f },   // CTERM
    { "CLYS", " NZ ", +0.6381f },   // CTERM
    { "CLYS", " O  ", -0.8252f },   // CTERM
    { "CLYS", " OXT", -0.8252f },   // CTERM
    { "CMET", " C  ", +0.8013f },   // CTERM
    { "CMET", " CA ", -0.1320f },   // CTERM
    { "CMET", " CB ", +0.0724f },   // CTERM
    { "CMET", " CE ", +0.1499f },   // CTERM
    { "CMET", " CG ", +0.1126f },   // CTERM
    { "CMET", " N  ", -0.1140f },   // CTERM
    { "CMET", " O  ", -0.8105f },   // CTERM
    { "CMET", " OXT", -0.8105f },   // CTERM
    { "CMET", " SD ", -0.2692f },   // CTERM
    { "CPHE", " C  ", +0.7660f },   // CTERM
    { "CPHE", " CA ", -0.0727f },   // CTERM
    { "CPHE", " CB ", -0.0073f },   // CTERM
    { "CPHE", " CD1", +0.0108f },   // CTERM
    { "CPHE", " CD2", +0.0108f },   // CTERM
    { "CPHE", " CE1", -0.0386f },   // CTERM
    { "CPHE", " CE2", -0.0386f },   // CTERM
    { "CPHE", " CG ", +0.0552f },   // CTERM
    { "CPHE", " CZ ", +0.0336f },   // CTERM
    { "CPHE", " N  ", -0.1140f },   // CTERM
    { "CPHE", " O  ", -0.8026f },   // CTERM
    { "CPHE", " OXT", -0.8026f },   // CTERM
    { "CPRO", " C  ", +0.6631f },   // CTERM
    { "CPRO", " CA ", -0.0560f },   // CTERM
    { "CPRO", " CB ", +0.0219f },   // CTERM
    { "CPRO", " CD ", +0.1096f },   // CTERM
    { "CPRO", " CG ", +0.0810f },   // CTERM
    { "CPRO", " N  ", -0.2802f },   // CTERM
    { "CPRO", " O  ", -0.7697f },   // CTERM
    { "CPRO", " OXT", -0.7697f },   // CTERM
    { "CSER", " C  ", +0.8113f },   // CTERM
    { "CSER", " CA ", -0.1418f },   // CTERM
    { "CSER", " CB ", +0.2749f },   // CTERM
    { "CSER", " N  ", -0.1140f },   // CTERM
    { "CSER", " O  ", -0.8132f },   // CTERM
    { "CSER", " OG ", -0.2040f },   // CTERM
    { "CSER", " OXT", -0.8132f },   // CTERM
    { "CTHR", " C  ", +0.7810f },   // CTERM
    { "CTHR", " CA ", -0.1213f },   // CTERM
    { "CTHR", " CB ", +0.3103f },   // CTERM
    { "CTHR", " CG2", -0.0095f },   // CTERM
    { "CTHR", " N  ", -0.1140f },   // CTERM
    { "CTHR", " O  ", -0.8044f },   // CTERM
    { "CTHR", " OG1", -0.2377f },   // CTERM
    { "CTHR", " OXT", -0.8044f },   // CTERM
    { "CTRP", " C  ", +0.7658f },   // CTERM
    { "CTRP", " CA ", -0.0812f },   // CTERM
    { "CTRP", " CB ", +0.0252f },   // CTERM
    { "CTRP", " CD1", +0.0235f },   // CTERM
    { "CTRP", " CD2", +0.1078f },   // CTERM
    { "CTRP", " CE2", +0.1222f },   // CTERM
    { "CTRP", " CE3", -0.0346f },   // CTERM
    { "CTRP", " CG ", -0.0796f },   // CTERM
    { "CTRP", " CH2", +0.0381f },   // CTERM
    { "CTRP", " CZ2", -0.1027f },   // CTERM
    { "CTRP", " CZ3", -0.0780f },   // CTERM
    { "CTRP", " N  ", -0.1140f },   // CTERM
    { "CTRP", " NE1", +0.0097f },   // CTERM
    { "CTRP", " O  ", -0.8011f },   // CTERM
    { "CTRP", " OXT", -0.8011f },   // CTERM
    { "CTYR", " C  ", +0.7817f },   // CTERM
    { "CTYR", " CA ", -0.0923f },   // CTERM
    { "CTYR", " CB ", +0.0228f },   // CTERM
    { "CTYR", " CD1", -0.0142f },   // CTERM
    { "CTYR", " CD2", -0.0142f },   // CTERM
    { "CTYR", " CE1", -0.0785f },   // CTERM
    { "CTYR", " CE2", -0.0785f },   // CTERM
    { "CTYR", " CG ", +0.0243f },   // CTERM
    { "CTYR", " CZ ", +0.3395f },   // CTERM
    { "CTYR", " N  ", -0.1140f },   // CTERM
    { "CTYR", " O  ", -0.8070f },   // CTERM
    { "CTYR", " OH ", -0.1626f },   // CTERM
    { "CTYR", " OXT", -0.8070f },   // CTERM
    { "CVAL", " C  ", +0.8350f },   // CTERM
    { "CVAL", " CA ", -0.2000f },   // CTERM
    { "CVAL", " CB ", +0.2248f },   // CTERM
    { "CVAL", " CG1", -0.0556f },   // CTERM
    { "CVAL", " CG2", -0.0556f },   // CTERM
    { "CVAL", " N  ", -0.1140f },   // CTERM
    { "CVAL", " O  ", -0.8173f },   // CTERM
    { "CVAL", " OXT", -0.8173f },   // CTERM
    { "HIS", " C  ", +0.5973f },   // STD
    { "HIS", " CA ", +0.0779f },   // STD
    { "HIS", " CB ", +0.0660f },   // STD
    { "HIS", " CD2", -0.0345f },   // STD
    { "HIS", " CE1", +0.3070f },   // STD
    { "HIS", " CG ", +0.1868f },   // STD
    { "HIS", " N  ", -0.1438f },   // STD
    { "HIS", " ND1", -0.5432f },   // STD
    { "HIS", " NE2", +0.0544f },   // STD
    { "HIS", " O  ", -0.5679f },   // STD
};
static constexpr int FF14SB_LUMPED_COUNT = 618;
} // namespace ff14sb_lumped
#endif // FF14SB_LUMPED_CHARGES_H
