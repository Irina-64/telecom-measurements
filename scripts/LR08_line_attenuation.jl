using Statistics

const NSAMP = 2000
const DT = 0.001
const F0 = 5.0
const A_IN = 2.0
const ALPHA_DB_PER_KM = 0.2
const LENGTH_KM = 50.0
const SEED = 101

A_theory_dB(alpha, L) = alpha * L

function amp_at_freq(x, f::Float64, fs::Float64, n::Int)
    k = round(Int, f * n / fs)
    X = 0.0 + 0.0im
    for m in 0:n-1
        X += x[m+1] * cispi(-2*k*m/n)
    end
    return 2 * abs(X) / n
end

function run_lr08(seed::Integer)
    engee.set_param!("LR08_line_attenuation/MeasNoise", "Seed" => string(seed))
    engee.run("LR08_line_attenuation")
    data = copy(measured_line_output.value)
    fs = 1 / DT
    a_out = amp_at_freq(data, F0, fs, NSAMP)
    return (data=data, a_out=a_out)
end

r1 = run_lr08(SEED)
r2 = run_lr08(SEED)
r3 = run_lr08(SEED + 1000)
r4 = run_lr08(SEED)

EQ_1_2 = r1.data == r2.data
EQ_1_4 = r1.data == r4.data
EQ_1_3 = r1.data == r3.data
EQ_2_4 = r2.data == r4.data

A_meas_dB = 20 * log10(A_IN / r1.a_out)
A_theory = A_theory_dB(ALPHA_DB_PER_KM, LENGTH_KM)

println("a_out=", r1.a_out, " A_meas_dB=", A_meas_dB, " A_theory_dB=", A_theory)
println("EQ_1_2=", EQ_1_2, " EQ_1_4=", EQ_1_4, " EQ_1_3=", EQ_1_3, " EQ_2_4=", EQ_2_4)

ratio = A_meas_dB / A_theory
PASS = EQ_1_2 && EQ_1_4 && !EQ_1_3 && EQ_2_4 && 0.5 <= ratio <= 2.0
println("PASS=", PASS)
