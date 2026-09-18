using Statistics

const NSAMP = 2000
const DT = 0.001
const TAU = 0.05
const TEST_FREQS = [1.0, 3.0, 10.0]
const A_IN = 1.0
const SEED = 101

Htheory(f, tau) = 1 / sqrt(1 + (2*pi*f*tau)^2)

function amp_at_freq(x, f::Float64, fs::Float64, n::Int)
    k = round(Int, f * n / fs)
    X = 0.0 + 0.0im
    for m in 0:n-1
        X += x[m+1] * cispi(-2*k*m/n)
    end
    return 2 * abs(X) / n
end

function run_lr07(seed::Integer)
    engee.set_param!("LR07_transfer_characteristic/MeasNoise", "Seed" => string(seed))
    engee.run("LR07_transfer_characteristic")
    data = copy(measured_output.value)
    fs = 1 / DT
    amps = [amp_at_freq(data, f, fs, NSAMP) for f in TEST_FREQS]
    return (data=data, amps=amps)
end

r1 = run_lr07(SEED)
r2 = run_lr07(SEED)
r3 = run_lr07(SEED + 1000)
r4 = run_lr07(SEED)

EQ_1_2 = r1.data == r2.data
EQ_1_4 = r1.data == r4.data
EQ_1_3 = r1.data == r3.data
EQ_2_4 = r2.data == r4.data

H_meas = r1.amps ./ A_IN
H_theory = [Htheory(f, TAU) for f in TEST_FREQS]

for i in 1:length(TEST_FREQS)
    println("f=", TEST_FREQS[i], " H_meas=", H_meas[i], " H_theory=", H_theory[i])
end
println("EQ_1_2=", EQ_1_2, " EQ_1_4=", EQ_1_4, " EQ_1_3=", EQ_1_3, " EQ_2_4=", EQ_2_4)

ratios = H_meas ./ H_theory
PASS = EQ_1_2 && EQ_1_4 && !EQ_1_3 && EQ_2_4 && all(0.5 .<= ratios .<= 2.0)
println("PASS=", PASS)
