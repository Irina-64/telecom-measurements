using Statistics
using FFTW

const MODEL_NAME = "LR05_spectral_windowing"
const RANDOM_BLOCK_PATH = "LR05_spectral_windowing/RandomNoise"
const BASE_SEED = 101
const ALT_SEED = BASE_SEED + 1000
const EXPECTED_LENGTH = 1000
const FS = 100.0
const PROBE_BIN = 70

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
    
    function window_metrics(x)
        N = length(x)
            half = 1:div(N,2)
                Yr = fft(x)
                    Pr = abs2.(Yr) ./ N
    k0r = argmax(Pr[half .+ 1])
    f_rect = k0r * FS / N
    w = [0.5*(1-cos(2*pi*(n-1)/(N-1))) for n in 1:N]
    xw = x .* w
    Yw = fft(xw)
    cg = sum(w)/N
        Pw = abs2.(Yw) ./ (N * cg^2 * N)
            k0w = argmax(Pw[half .+ 1])
    f_hann = k0w * FS / N
        leakage_rect = Pr[PROBE_BIN]
    leakage_hann = Pw[PROBE_BIN]
        leakage_ratio = leakage_hann / leakage_rect
            return (n=N, f_rect=f_rect, f_hann=f_hann, leakage_rect=leakage_rect, leakage_hann=leakage_hann, leakage_ratio=leakage_ratio)
            end
            
            function check_numeric_acceptance(m)
                @assert isapprox(m.f_rect, 5.25; atol=0.15) "Частота (прямоугольное окно) вне допуска"
                    @assert isapprox(m.f_hann, 5.25; atol=0.15) "Частота (Hann) вне допуска"
                        @assert m.leakage_ratio < 0.01 "Окно Hann не дало ожидаемого снижения утечки"
                            return true
end

println("MODEL = ", MODEL_NAME)
println("SEEDS = ", (BASE_SEED, BASE_SEED, ALT_SEED, BASE_SEED))

set_seed!(BASE_SEED)
x1 = run_model_and_extract()
m1 = window_metrics(x1)

set_seed!(BASE_SEED)
x2 = run_model_and_extract()
m2 = window_metrics(x2)

set_seed!(ALT_SEED)
x3 = run_model_and_extract()
m3 = window_metrics(x3)

set_seed!(BASE_SEED)
x4 = run_model_and_extract()
m4 = window_metrics(x4)

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
