#
# The one hand-curated ingredient `tools/buildCapabilityRegister.jl` reads before writing CAPABILITIES.md.
# Everything else in the register is re-derived mechanically from examples/examples.jl and the `Last
# successful` dates already recorded in examples/*.jl; this file holds only what cannot be derived --
# a judgement call, an abbreviation, or a title examples.jl itself never prints.
#
# Growing this file is always safe and never urgent: a branch with no entry here still gets a correct,
# honest row (just without a restriction hint, a keyword synonym, or an explicit module name) -- nothing
# here is required for the register to be correct, only for it to be a little more helpful.
#

# Section titles for branch-letter groups that examples/examples.jl itself never prints a header for --
# 'Q' (muonic atoms) is folded into the 'P' (photon-scattering) menu entry in that file, and 'S', 'T' have
# no header println at all.
const GROUP_TITLES = Dict(
    "Q" => "Muonic atoms",
    "S" => "Stark / line-broadening (semiempirical) computations",
    "T" => "Three-body recombination",
)

# Module overrides, applied AFTER automatic extraction and only where auto-extraction cannot find a clean
# "the X module" / "(module X)" phrase, verified against src/ rather than guessed.
const MODULE_OVERRIDE = Dict(
    "Ai" => "AtomicState.RasExpansion",
    "Ao" => "Basics (AverageLevelRoute / RotationRoute -- the AL/EOL SCF routes)",
    "Ar" => "Hamiltonian.performCIwithSecondOrderQ",
    "Ba" => "Basics (dipole / multipole amplitude routines)",
)

# Restriction hints -- one line, deliberately not detailed (see CAPABILITIES.md's own preamble). Keyed by
# branch code. A missing key means no restriction is currently flagged, not that none exists.
const RESTRICTIONS = Dict(
    "Ai" => "cost grows fast with layer count and nuclear charge -- see `docs/` or ask before a 4th layer on a heavy ion.",
    "Cb" => "a subshell-factor convention fix landed in 2026 (`CLAUDE.md` Rule 18) -- a number computed before that fix is not comparable to one after it.",
    "De" => "a modelling edge case for this process is under active review; check for a more recent date before relying on a result close to a known difficult case (e.g. near-degenerate continua).",
    "Df" => "this file carries several branches that are still undated -- check the date of the SPECIFIC branch you copy, not just the file.",
    "Dl" => "at least one tested case has come in well below the high-energy Bethe benchmark; treat an absolute rate as a first estimate until checked against an independent source.",
    "Ja" => "the free-free contribution near zero frequency is an open point; do not rely on this model there without checking for a newer date.",
    "Id" => "explicitly marked [POSTPONED until 2027] in examples.jl itself -- do not use for a real calculation.",
    "Ic" => "shares the StrongField module with Id; check this branch's own date independently before trusting it.",
    "Ib" => "shares the StrongField module with Id; check this branch's own date independently before trusting it.",
)

# Keyword synonyms -- only for abbreviations or alternate names NOT already a substring of the entry's own
# description, so the plain-text search tool still finds the entry under the term a user is likely to type.
const SYNONYMS = Dict(
    "Ai" => ["RAS", "restricted active space", "correlation layer"],
    "Cb" => ["HFS", "hyperfine", "A constant", "B constant"],
    "Cr" => ["HFS", "Breit-Rabi", "Paschen-Back"],
    "Do" => ["HFS"],
    "De" => ["Auger"],
    "Df" => ["DR"],
    "Dx" => ["DR", "RR"],
    "Dy" => ["Fano"],
    "Fa" => ["cascade", "Auger cascade"],
    "Fd" => ["cascade"],
    "Ff" => ["cascade", "KLL"],
    "Fm" => ["cascade", "KLL"],
    "Ja" => ["average atom", "plasma", "finite temperature"],
    "Jc" => ["Saha-Boltzmann", "LTE", "ionization balance"],
    "Hc" => ["Zeff", "Rydberg", "quick estimate"],
    "Ma" => ["ForPedestrians", "quick start", "beginner"],
    "Mb" => ["ForPedestrians", "quick start", "beginner"],
    "Mc" => ["ForPedestrians", "quick start", "beginner"],
    "Md" => ["ForPedestrians", "quick start", "beginner"],
    "Me" => ["ForPedestrians", "quick start", "beginner"],
    "Mf" => ["ForPedestrians", "quick start", "beginner"],
    "Mg" => ["ForPedestrians", "quick start", "beginner"],
    "Mh" => ["ForPedestrians", "quick start", "beginner"],
    "Mi" => ["ForPedestrians", "quick start", "beginner"],
    "Id" => ["ATI", "strong field", "SFA"],
    "Ic" => ["ATI", "strong field", "SFA"],
    "Ib" => ["ATI", "strong field", "SFA"],
)
