using Statistics

const NSAMP = 2000
const BASELINE = 20.0
const DELTA_L = 15.0
const SEED = 101

function run_lr13(seed::Integer)
    engee.set_param!("LR13_metrics_logs_traces/Jitter", "Seed" => string(seed))
    engee.run("LR13_metrics_logs_traces")
    metric = collect(Float64, metric_latency.value)
    logsig = collect(Float64, log_fault_signal.value)
    fault_flag = logsig .> 0.5
    baseline_mean = mean(metric[.!fault_flag])
    fault_mean = mean(metric[fault_flag])
    r = cor(metric, Float64.(fault_flag))
    return (metric=metric, logsig=logsig, baseline_mean=baseline_mean, fault_mean=fault_mean, r=r)
end

r1 = run_lr13(SEED)
r2 = run_lr13(SEED)
r3 = run_lr13(SEED + 1000)
r4 = run_lr13(SEED)

EQ_m_1_2 = r1.metric == r2.metric
EQ_m_1_4 = r1.metric == r4.metric
EQ_m_1_3 = r1.metric == r3.metric
EQ_m_2_4 = r2.metric == r4.metric

EQ_l_1_2 = r1.logsig == r2.logsig
EQ_l_1_4 = r1.logsig == r4.logsig
EQ_l_1_3 = r1.logsig == r3.logsig
EQ_l_2_4 = r2.logsig == r4.logsig

diff_means = r1.fault_mean - r1.baseline_mean

println("baseline_mean=", r1.baseline_mean, " fault_mean=", r1.fault_mean, " diff=", diff_means, " DELTA_L=", DELTA_L)
println("correlation r=", r1.r)
println("EQ_metric: ", EQ_m_1_2, " ", EQ_m_1_4, " ", EQ_m_1_3, " ", EQ_m_2_4)
println("EQ_log: ", EQ_l_1_2, " ", EQ_l_1_4, " ", EQ_l_1_3, " ", EQ_l_2_4)

ratio = diff_means / DELTA_L
PASS = EQ_m_1_2 && EQ_m_1_4 && !EQ_m_1_3 && EQ_m_2_4 &&
       EQ_l_1_2 && EQ_l_1_4 && EQ_l_1_3 && EQ_l_2_4 &&
       0.8 <= ratio <= 1.2 && r1.r > 0.7
println("PASS=", PASS)
