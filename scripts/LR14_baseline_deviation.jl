using Statistics

const NSAMP = 2000
const TRAIN_N = 1000
const ANOMALY_START = 10.0
const DT = 0.01
const K_SIGMA = 3.0
const SEED = 101

function run_lr14(seed::Integer)
    engee.set_param!("LR14_baseline_deviation/Jitter", "Seed" => string(seed))
    engee.run("LR14_baseline_deviation")
    data = collect(Float64, metric_value.value)
    train = data[1:TRAIN_N]
    baseline_mean = mean(train)
    baseline_std = std(train)
    threshold = baseline_mean + K_SIGMA * baseline_std
    flag = data .> threshold
    is_anomaly_period = [((k-1)*DT) >= ANOMALY_START for k in 1:NSAMP]
    tp = count(flag .& is_anomaly_period)
    fp = count(flag .& .!is_anomaly_period)
    n_anomaly = count(is_anomaly_period)
    n_normal = NSAMP - n_anomaly
    tpr = tp / n_anomaly
    fpr = fp / n_normal
    return (data=data, baseline_mean=baseline_mean, baseline_std=baseline_std, threshold=threshold, tpr=tpr, fpr=fpr)
end

r1 = run_lr14(SEED)
r2 = run_lr14(SEED)
r3 = run_lr14(SEED + 1000)
r4 = run_lr14(SEED)

EQ_1_2 = r1.data == r2.data
EQ_1_4 = r1.data == r4.data
EQ_1_3 = r1.data == r3.data
EQ_2_4 = r2.data == r4.data

println("baseline_mean=", r1.baseline_mean, " baseline_std=", r1.baseline_std, " threshold=", r1.threshold)
println("TPR (detection rate)=", r1.tpr, " FPR (false positive rate)=", r1.fpr)
println("EQ_1_2=", EQ_1_2, " EQ_1_4=", EQ_1_4, " EQ_1_3=", EQ_1_3, " EQ_2_4=", EQ_2_4)

PASS = EQ_1_2 && EQ_1_4 && !EQ_1_3 && EQ_2_4 && r1.tpr > 0.9 && r1.fpr < 0.05
println("PASS=", PASS)
