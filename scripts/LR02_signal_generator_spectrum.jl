using Statistics
using FFTW

const MODEL_NAME = "LR02_signal_generator_spectrum"
const RANDOM_BLOCK_PATH = "LR02_signal_generator_spectrum/RandomNoise"
const BASE_SEED = 101
const ALT_SEED = BASE_SEED + 1000
const EXPECTED_LENGTH = 1000

function run_model_and_extract()
    engee.run(MODEL_NAME)
        x = Float64.(collect(generated_signal).value)
            @assert !isempty(x) "Выход generated_signal пуст"
                @assert length(x) == EXPECTED_LENGTH "Неожиданная длина выхода"
                    @assert all(isfinite, x) "Выход содержит NaN или Inf"
                        return x
                        end

function set_seed!(seed::Integer)
    engee.set_param!(RANDOM_BLOCK_PATH, "Seed" => string(seed))
    end

function spectral_metrics(x)
    N = length(x)
    fs = 100.0
        Y = fft(x)
            P = abs2.(Y) ./ N
    half = 1:div(N,2)
        k0 = argmax(P[half .+ 1])
            f0_hat = k0 * fs / N
    A0_hat = 2*abs(Y[k0+1])/N
        P_total = mean(x.^2)
            P_signal = A0_hat^2/2
                P_noise = P_total - P_signal
                    SNR_dB = 10*log10(P_signal/P_noise)
                        return (n=N, f0=f0_hat, A0=A0_hat, SNR=SNR_dB)
                        end

function check_numeric_acceptance(m)
    @assert isapprox(m.f0, 5.0; atol=0.05) "Частота вне допуска"
        @assert 1.8 <= m.A0 <= 2.2 "Амплитуда вне допуска"
            @assert 14.0 <= m.SNR <= 20.0 "SNR вне допуска"
                return true
end

println("MODEL = ", MODEL_NAME)
println("SEEDS = ", (BASE_SEED, BASE_SEED, ALT_SEED, BASE_SEED))

set_seed!(BASE_SEED)
x1 = run_model_and_extract()
m1 = spectral_metrics(x1)

set_seed!(BASE_SEED)
x2 = run_model_and_extract()
m2 = spectral_metrics(x2)

set_seed!(ALT_SEED)
x3 = run_model_and_extract()
m3 = spectral_metrics(x3)

set_seed!(BASE_SEED)
x4 = run_model_and_extract()
m4 = spectral_metrics(x4)

eq_1_2 = x1 == x2
eq_1_4 = x1 == x4
eq_1_3 = x1 == x3
eq_2_4 = x2 == x4

@assert eq_1_2 "Запуски 1 и 2 при одинаковом seed различаются"
@assert eq_1_4 "После восстановления seed выход не воспроизведён"
@assert !eq_1_3 "Изменение seed не изменило полный выходной массив"
@assert eq_2_4 "Запуски 2 и 4 при одинаковом seed различаются"

check_numeric_acceptance(m1)
check_numeric_acceptance(m2)
check_numeric_acceptance(m4)

set_seed!(BASE_SEED)
engee.save(MODEL_NAME, MODEL_NAME * ".engee"; force=true)

println("RUN_1 = ", m1)
println("RUN_2 = ", m2)
println("RUN_3 = ", m3)
println("RUN_4 = ", m4)
println("EQ_1_2 = ", eq_1_2)
println("EQ_1_4 = ", eq_1_4)
println("EQ_1_3 = ", eq_1_3)
println("EQ_2_4 = ", eq_2_4)
println("PASS: модель и критерии варианта 1 проверены; параметры восстановлены")
