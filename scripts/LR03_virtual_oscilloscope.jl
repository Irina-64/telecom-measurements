using Statistics

const MODEL_NAME = "LR03_virtual_oscilloscope"
const RANDOM_BLOCK_PATH = "LR03_virtual_oscilloscope/RandomNoise"
const BASE_SEED = 101
const ALT_SEED = BASE_SEED + 1000
const EXPECTED_LENGTH = 2000
const THRESHOLD = 1.0

function run_model_and_extract()
    engee.run(MODEL_NAME)
        df = collect(oscilloscope_trace)
            t = Float64.(df.time)
                x = Float64.(df.value)
                    @assert !isempty(x) "Выход oscilloscope_trace пуст"
    @assert length(x) == EXPECTED_LENGTH "Неожиданная длина выхода"
    @assert all(isfinite, x) "Выход содержит NaN или Inf"
    return t, x
end

function set_seed!(seed::Integer)
    engee.set_param!(RANDOM_BLOCK_PATH, "Seed" => string(seed))
    end

function rising_edges(t, x, thr)
    ts = Float64[]
        for i in 2:length(x)
                if x[i-1] < thr && x[i] >= thr
            frac = (thr - x[i-1]) / (x[i] - x[i-1])
                        push!(ts, t[i-1] + frac*(t[i]-t[i-1]))
                                end
    end
        return ts
end

function quantized_edges(t, x, thr)
    ts = Float64[]
        for i in 2:length(x)
                if x[i-1] < thr && x[i] >= thr
                            push!(ts, t[i])
                                    end
                                        end
                                            return ts
                                            end

function oscilloscope_metrics(t, x)
    edges = rising_edges(t, x, THRESHOLD)
        edges_q = quantized_edges(t, x, THRESHOLD)
            periods = diff(edges)
    T_hat = mean(periods)
        f_hat = 1/T_hat
    jitter_std = std(periods)
        ts_err = abs.(edges_q .- edges)
    Vpp = maximum(x) - minimum(x)
        high_level = mean(x[x .> THRESHOLD])
            low_level = mean(x[x .<= THRESHOLD])
                return (n=length(x), nedges=length(edges), T=T_hat, f=f_hat, jitter=jitter_std, ts_err_mean=mean(ts_err), ts_err_max=maximum(ts_err), Vpp=Vpp, high=high_level, low=low_level)
                end

function check_numeric_acceptance(m)
    @assert 8 <= m.nedges <= 10 "Число фронтов вне диапазона"
        @assert isapprox(m.T, 0.2; atol=0.01) "Период вне допуска"
            @assert m.jitter < 0.0005 "Джиттер превышает допустимый"
                @assert m.ts_err_max <= 0.001 "Ошибка временной метки превышает шаг дискретизации"
                    @assert 1.8 <= m.high <= 2.2 "Уровень логической единицы вне допуска"
    return true
end

println("MODEL = ", MODEL_NAME)
println("SEEDS = ", (BASE_SEED, BASE_SEED, ALT_SEED, BASE_SEED))

set_seed!(BASE_SEED)
t1, x1 = run_model_and_extract()
m1 = oscilloscope_metrics(t1, x1)

set_seed!(BASE_SEED)
t2, x2 = run_model_and_extract()
m2 = oscilloscope_metrics(t2, x2)

set_seed!(ALT_SEED)
t3, x3 = run_model_and_extract()
m3 = oscilloscope_metrics(t3, x3)

set_seed!(BASE_SEED)
t4, x4 = run_model_and_extract()
m4 = oscilloscope_metrics(t4, x4)

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
