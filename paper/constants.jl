## General constants 
const EXTENSION = "pdf"
const NIST = true
const FIGURES = "paper/figures"
const ALPHABET = join('a':'z')

# The constants below are obtained in the MATLAB script "burger_aaa.m".
const U_AAA_1 = parse.(ComplexF64, readdlm("paper/data/U_AAA_1.dat"))
const U_AAA_2 = parse.(ComplexF64, readdlm("paper/data/U_AAA_2.dat"))
const U_AAA_3 = parse.(ComplexF64, readdlm("paper/data/U_AAA_3.dat"))
const U_AAA_4 = parse.(ComplexF64, readdlm("paper/data/U_AAA_4.dat"))
const U_AAA = [U_AAA_1, U_AAA_2, U_AAA_3, U_AAA_4]
const X_AAA = vec(readdlm("paper/data/X_AAA.dat"))
const Y_AAA = vec(readdlm("paper/data/Y_AAA.dat"))
const T_AAA = vec(readdlm("paper/data/T_AAA.dat"))
const μ_AAA = readdlm("paper/data/mu_AAA.dat")[1]

# The constants below are obtained in the MATLAB script "burger_tracking_poles.m".
const T_AAA_POLES = vec(readdlm("paper/data/T_AAA_POLES.dat"))
const AAA_POLES = parse.(ComplexF64, readdlm("paper/data/AAA_POLES.dat"))