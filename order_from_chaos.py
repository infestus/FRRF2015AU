#! /usr/bin/python3
from gi import pygtkcompat
import numpy
import math
import math_plot as mp
import sys
import seqCount as sq
pygtkcompat.enable() 
pygtkcompat.enable_gtk(version='3.0')
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from multiprocessing import Process

def log_time_diff(d, big = 2000, low = 400, b=5):
    log_inc = ((mp.log(big,b)-mp.log(low,b))/mp.log(low, b))*100
    big_scaled = d[5][big]/mp.rna_running_time(big,b) 
    low_scaled = d[5][low]/mp.rna_running_time(low,b)
    time_inc = ((big_scaled-low_scaled)/low_scaled)*100
    print("Percentage logarithmic increase: " + str(log_inc) + " vs " + "percentage time increase: " + str(time_inc))
    print("Scaled percentage difference: " + str(low_scaled/big_scaled))

def log_incperc(big, low, b = 5):
    return ((math.log(big,b)-math.log(low,b))/math.log(low, b))*100


def data_incperc(d, big, low, b = 5):
    big_scaled = d[5][big]/rna_running_time(big,b) 
    low_scaled = d[5][low]/rna_running_time(low,b)
    return ((big_scaled-low_scaled)/low_scaled)*100

def median(l):
    return l[int(len(l)/2)-1]

def rna_running_time(n, b = 5):
    return math.pow(n, 3)/(round(math.log(n,b)+ 0.000001))

def create_burn_data_dict(filename, burn):
    file = open(filename, "r")
    d = {}
    for line in file:
        l = line.split(";")
        n = int(l[0])
        q = int(l[1])
        burn_in = (1 == int(l[2].split(":")[1]))
        t = float(l[3])
        if (burn == burn_in):
            if q not in d:
                d[q] = {n : [t]}
            else:
                if (n in d[q]):
                    d[q][n].append(t)
                else:
                    d[q][n] = [t]
    return d

def create_data_points_dict(filename):
    file = open(filename, "r")
    file.readline()
    d = {}
    for line in file: 
        l = line.split(";")
        n = int(l[0])
        q = int(l[1])
        t = float(l[2])
        if q not in d:
            d[q] = {n : [t]}
        else:
            if n in d[q]: 
                d[q][n].append(t)
            else: 
                d[q][n] = [t]
    return d

def create_q_indexed_dict(filename):
    data = create_data_points_dict(filename)
    d = {}
    for q in data:
        for n in data[q]:
            s = math.fsum(data[q][n])
            mean = s/len(data[q][n])
            if q not in d:
                d[q] = {n : mean}
            else: 
                d[q][n] = mean
    return d;

def create_q_indexed_dict_median(filename):
    data = create_data_points_dict(filename)
    d = {}
    for q in data:
        for n in data[q]:
            sdata = sorted(data[q][n])
            if q not in d:
                d[q] = {n : (sdata[int(len(sdata)/2)-1])}
            else: 
                d[q][n] = sdata[int(len(sdata)/2)-1]
    return d

def create_n_indexed_dict(filename):
    data = create_data_points_dict(filename)
    d = {}
    for q in data:
        for n in data[q]:
            s = math.fsum(data[q][n])
            mean = s/len(data[q][n])
            if n not in d:
                d[n] = {q : mean}
            else: 
                d[n][q] = mean
    return d;

def create_n_indexed_dict_median(filename):
    data = create_data_points_dict(filename)
    d = {}
    for q in data:
        for n in data[q]:
            sdata = sorted(data[q][n])
            if n not in d:
                d[n] = {q : (sdata[int(len(sdata)/2)-1])}
            else: 
                d[n][q] = sdata[int(len(sdata)/2)-1]
    return d;

main_data_file = "./data/new_2k.dat"

data_points = create_data_points_dict(main_data_file)

data_dict_q = create_q_indexed_dict_median(main_data_file)
data_dict = create_n_indexed_dict_median(main_data_file)
q_scaled_dict = create_n_indexed_dict_median("./data/1k_unique.dat")
time_scaled_dict = create_q_indexed_dict_median("./data/new_6k.dat")
big_6k_dict = create_q_indexed_dict_median("./data/new_6k.dat")
data_dict_q_35h = create_q_indexed_dict_median("./data/new_35h.dat")
uneven_2k = create_q_indexed_dict_median("./data/uneven2k.dat")
burn_ins_dict = create_burn_data_dict("./data/burn_test_3.dat", True)
no_burn_ins_dict = create_burn_data_dict("./data/burn_test_3.dat", False)
dynamic_dict = create_q_indexed_dict_median("./data/dynamic_unique6k.dat")
static_dict = create_q_indexed_dict_median("./data/static_unique6k.dat")

plt.xticks(rotation=270)
plt.tight_layout()

fig1 = plt.figure(1)
lin_plt = fig1.add_subplot(111)
lin_plt.grid(True)
lin_plt.set_xlabel("Input Size")
lin_plt.set_ylabel("Time in seconds") 
for q in data_dict_q:
    n_s = []
    v_s = []
    for n in sorted(data_dict_q[q]):
        n_s.append(n)
        v_s.append(data_dict_q[q][n])
    if q==0:
        lin_plt.plot(n_s, v_s, '-.', label=("Nussunov"), linewidth=2.0)
    elif q==-1:
        # lin_plt.plot(n_s, v_s, 'go', label=("Q: Dynamic"))
        continue; 
    else:
        lin_plt.plot(n_s, v_s, label=("Q:" + str(q)))
lin_plt.legend(bbox_to_anchor=(0.5, 1.1), loc='upper center', borderaxespad=0., ncol=3, fancybox=True, shadow=True)
fig1.savefig("./plots/linear.png", bbox_inches='tight', pad_inches=0.2)


fig2 = plt.figure(2)
rel_plt = fig2.add_subplot(111) 
rel_plt.grid(True)
rel_plt.set_xlabel("Input Size")
rel_plt.set_ylabel("Relative Performance vs Nussinov\n(higher is better)")
for q in data_dict_q:
    n_s = []
    v_s = []
    if q==0:
        continue
    else:
        for n in sorted(data_dict_q[q]):
            if n > 1000:
                break
            n_s.append(n)
            v_s.append(data_dict_q[0][n]/data_dict_q[q][n])
        rel_plt.plot(n_s, v_s, '-', label=("Q:" + str(q)))
for b in dynamic_dict:
    n_s = []
    v_s = []
    if b==0:
        continue
    else:
        for n in sorted(dynamic_dict[b]):
            if (n > 1000):
                break
            n_s.append(n)
            v_s.append(dynamic_dict[0][n]/dynamic_dict[b][n])
        rel_plt.plot(n_s, v_s, '-.', label=("b:" + str(b)))

rel_plt.legend(bbox_to_anchor=(0.5, 1.3), loc='upper center', borderaxespad=0., ncol=3, fancybox=True, shadow=True)
fig2.savefig("./plots/relative.png", bbox_inches='tight', pad_inches=0.2)


fig3 = plt.figure(3)
dots = fig3.add_subplot(111)
dots.grid(True)
dots.set_xlabel("Input Size")
dots.set_ylabel("Time in seconds")
for q in data_points:
    n_s = []
    v_s = []
    if q==5 or q==0:
        for n in sorted(data_points[q]):
            n_s = n_s + [n for a in range(len(data_points[q][n]))]
            v_s = v_s + data_points[q][n]
    if (q==5):
        dots.plot(n_s, v_s, 'x', label=("Q:5"))
    elif (q==0): 
        dots.plot(n_s, v_s, '.', label=("Nussinov"))
    else:
        continue
dots.legend(bbox_to_anchor=(0.5, 1.1), loc='upper center', borderaxespad=0., ncol=3, fancybox=True, shadow=True)
sq = 0
mean = 0
stdlist_rna = {}
stdlist_nussinov = {}
for n in data_points[5]:

    ls = sorted(data_points[5][n])
    mean_rna = sum(ls)/len(ls)
    variances = [math.pow(i-mean_rna,2) for i in ls]
    variance_rna = sum(variances)/len(variances)
    rna_rstddiv = round((math.sqrt(variance_rna)/mean_rna)*100, 2)
    stdlist_rna[n] = rna_rstddiv

    ls = sorted(data_points[0][n])
    mean_nussinov = sum(ls)/len(ls)
    variances = [math.pow(i-mean_nussinov,2) for i in ls]
    variance_nussinov = sum(variances)/len(variances)
    nussinov_rstddiv = round((math.sqrt(variance_nussinov)/mean_nussinov)*100, 2)
    stdlist_nussinov[n] = nussinov_rstddiv
max_val = max(max(stdlist_rna, key=stdlist_rna.get), max(stdlist_nussinov, key=stdlist_nussinov.get))
max_val = max(stdlist_rna[max_val], stdlist_nussinov[max_val])
mean_stddeviation = 0
for n in stdlist_rna:
    mean_stddeviation += stdlist_rna[n]
    mean_stddeviation += stdlist_nussinov[n]
mean_stddeviation = round(mean_stddeviation/(2 * len(stdlist_rna)), 2)

tt = "Mean Std Deviation: " + str(mean_stddeviation) + "%\nMax Std Deviation:   " + str(max_val) + "%"
tt = tt.expandtabs()
fig3.text(0.365,0.85, tt, bbox=dict(facecolor="white", alpha=1, pad=20))

fig3.savefig("./plots/dots.png", bbox_inches='tight', pad_inches=0.2)


fig4 = plt.figure(4)
q_based = fig4.add_subplot(111)
q_based.grid(True)
q_based.set_xlabel("q")
q_based.set_ylabel("Time in seconds")
x_s = []
v_s = []
for q in q_scaled_dict[1000]:
    if q!=0:
        x_s.append(q)
        v_s.append(q_scaled_dict[1000][q])
q_based.plot(x_s, v_s, '-', label="n=1000")
q_based.legend(bbox_to_anchor=(0.5, 1.1), loc='upper center', borderaxespad=0., ncol=3, fancybox=True, shadow=True)
q_based.set_ylim([0, 50])
fig4.savefig("./plots/q_scaling.png", bbox_inches='tight', pad_inches=0.2)


# time_scaled_dict = create_q_indexed_dict_median("./data/all_unique.dat")

fig5 = plt.figure(5)
time_scaled = fig5.add_subplot(111)
time_scaled.grid(True)
time_scaled.set_xlabel("Input Size")
time_scaled.set_ylabel("Normalized running times")
n_s = []
static_q_values_slow = []
static_q_values_fast = []
nussinov_values = []
dynamic_q_values_slow = []
dynamic_q_values_fast = []
b = 5
q_val = 5
min_x = sorted(time_scaled_dict[q_val])[0]
max_x = sorted(time_scaled_dict[q_val])[-1]
for n in sorted(time_scaled_dict[q_val]):
    static_slow_scaled = (time_scaled_dict[q_val][n])/(math.pow(n,3))
    dyn_slow_scaled = (dynamic_dict[b][n])/math.pow(n,3)
    static_fast_scaled = (dynamic_dict[q_val][n])/(math.pow(n,3)/q_val)
    dyn_fast_scaled = (dynamic_dict[b][n])/rna_running_time(n, b)
    nus_scaled = (dynamic_dict[0][n])/(math.pow(n,3))
    static_q_values_slow.append(static_slow_scaled)
    static_q_values_fast.append(static_fast_scaled)
    dynamic_q_values_fast.append(dyn_fast_scaled)
    dynamic_q_values_slow.append(dyn_slow_scaled)
    nussinov_values.append(nus_scaled)
    n_s.append(n)
time_scaled.plot(n_s, nussinov_values, "b-.", label="Nussinov: O("+r'$n^3$' + ")-normalized", lw=2)
time_scaled.plot(n_s, static_q_values_slow, "g-.", label="q=5: O("+r'$n^3$' + ")-normalized", lw=3)
time_scaled.plot(n_s, static_q_values_fast, "g-", label="q=5: O("+r'$\frac{n^3}{5}$' + ")-normalized", lw=3)
time_scaled.plot(n_s, dynamic_q_values_slow, "ro-.", label="Dynamic: O("+r'$n^3$' + ")-normalized")
time_scaled.plot(n_s, dynamic_q_values_fast, "ro-", label="Dynamic: O("+r'$\frac{n^3}{log_5(n)}$' + ")-normalized")
time_scaled.legend(bbox_to_anchor=(0.5, 1.31), loc='upper center', borderaxespad=0., ncol=2, fancybox=True, shadow=True)
time_scaled.set_xlim([min_x,max_x])
fig5.savefig("./plots/time_scaled.png", bbox_inches='tight', pad_inches=0.2)


fig6 = plt.figure(6)
log_plot = fig6.add_subplot(111)
x_s = []
times_s = []
logs_s = []
inc_perc_dict = big_6k_dict
q_value = -1
min_n = sorted(inc_perc_dict[q_value])[0]
for n in sorted(inc_perc_dict[q_value]):
    if (n <= min_n):
        continue
    x_s.append(n)
    times_s.append(data_incperc(inc_perc_dict, n, min_n))
    logs_s.append(log_incperc(n, min_n))
log_plot.plot(x_s, times_s, label = "Four Russianss")
log_plot.plot(x_s, logs_s, label = "Log_5")
log_plot.legend(bbox_to_anchor=(0.5, 1.1), loc='upper center', borderaxespad=0., ncol=1, fancybox=True, shadow=True)
fig6.savefig("./plots/log_times.png", bbox_inches = 'tight', pad_inches = 0.2)


if (len(sys.argv) > 1):
    fig7 = plt.figure(7)
    h1 = fig7.add_subplot(111)
    h1.grid(True)
    seq_lens = sq.sequenceLengths("./sequences/rnacentral.fasta")
    cutoff = 0.05
    seq_lens = seq_lens[int((len(seq_lens)*cutoff)):int(len(seq_lens)-(len(seq_lens)*cutoff))]
    top = seq_lens[-1]
    bot = seq_lens[0]
    n, bins, patches = h1.hist(seq_lens, bins=[n for n in range(bot, top+1, 100)], range=[bot, top+1], facecolor='green', alpha=0.8)
    h1.set_xlim(bot, top)
    h1.set_xticks([n for i, n in enumerate(bins) if i % 2 == 1])
    h1.set_ylabel("Number of sequences")
    h1.set_xlabel("Sequence Lengths")
    s = sum(seq_lens)
    mean = int(s/len(seq_lens))
    median = seq_lens[int((len(seq_lens)/2)-1)]
    tt = "Mean:\t " + str(mean) + " nucleotides\nMedian:\t" + str(median) + " nucleotides"
    tt = tt.expandtabs()
    fig7.text(0.38,0.90, tt, bbox=dict(facecolor="white", alpha=1, pad=20))
    plt.xticks(rotation=270)
    #-------------------------
    lin = h1.twinx()
    x_s = [bot]
    y_s = [0]
    acc = 0
    for p in patches:
        acc = acc + p.get_height()
        x_s.append(p.get_x() + p.get_width())
        y_s.append((acc/len(seq_lens))*100)
    lin.plot(x_s, y_s, "-")
    lin.set_yticks([n for n in range(0, 101, 10)])
    lin.set_ylabel("Percentage of accumulated sequences")
    # lin.set_yticks([n for n in range(0, 100, 10)])
    fig7.savefig("./plots/lengthhistogram.png", bbox_inches='tight', pad_inches=0.2)


fig8 = plt.figure(8)
big = fig8.add_subplot(111)
big.grid(True)
y_s = []
x_s = []
for n in sorted(big_6k_dict[5]):
    y_s.append(big_6k_dict[5][n])
    x_s.append(n)
big.plot(x_s, y_s, "-", label="q = 5")
y_s = []
x_s = []
for n in sorted(big_6k_dict[0]):
    y_s.append(big_6k_dict[0][n])
    x_s.append(n)
big.plot(x_s, y_s, "-", label="Nussinov")
big.legend(bbox_to_anchor=(0.5, 1.1), loc='upper center', borderaxespad=0., ncol=1, fancybox=True, shadow=True)
fig8.savefig("./plots/biglinear.png", bbox_inches='tight', pad_inches=0.2)


fig9 = plt.figure(9)
uneven = fig9.add_subplot(111)
uneven.grid(True)
uneven.set_xlabel("Input Size")
uneven.set_ylabel("Time in seconds") 
for q in uneven_2k:
    n_s = []
    v_s = []
    for n in sorted(uneven_2k[q]):
        n_s.append(n)
        v_s.append(uneven_2k[q][n])
    if q==0:
        uneven.plot(n_s, v_s, '-.', label=("Nussunov"), linewidth=2.0)
    else:
        uneven.plot(n_s, v_s, label=("Q:" + str(q)))
uneven.legend(bbox_to_anchor=(0.5, 1.1), loc='upper center', borderaxespad=0., ncol=3, fancybox=True, shadow=True)
# fig9.savefig("./plots/uneven_linear.png", bbox_inches='tight', pad_inches=0.2)

fig10 = plt.figure(10)
burn_ins = fig10.add_subplot(111)
burn_ins.grid(True)
burn_ins.set_xlabel("Input Size")
burn_ins.set_ylabel("Burn-in quotient") 

for q in burn_ins_dict: 
    n_s = []
    v_s = []
    for n in sorted(burn_ins_dict[q]):
        n_s.append(n)
        v_s.append(median(burn_ins_dict[q][n])/median(no_burn_ins_dict[q][n]))
    # if (q == -1):
    #     burn_ins.plot(n_s, v_s, label="Dynamic")
    if (q == 0):
        burn_ins.plot(n_s, v_s, label="Nussinov")
    if (q > 0):
        burn_ins.plot(n_s, v_s, label="q " + str(q))
burn_ins.set_xlim([300,1000])
burn_ins.set_ylim([0.5,1.5])
burn_ins.legend(bbox_to_anchor=(0.5, 1.1), loc='upper center', borderaxespad=0., ncol=3, fancybox=True, shadow=True)
fig10.savefig("./plots/burn_ins.png", bbox_inches='tight', pad_inches=0.2)

fig11 = plt.figure(11)
dynamic_time = fig11.add_subplot(111)
dynamic_time.grid(True)
dynamic_time.set_xlabel("Input Size")
dynamic_time.set_ylabel("Time in seconds") 
for b in dynamic_dict:
    n_s = []
    v_s = []
    for n in sorted(dynamic_dict[b]):
        n_s.append(n)
        v_s.append(dynamic_dict[b][n])
    if (b == 0):
        dynamic_time.plot(n_s, v_s, "-.", label="Nussinov")
    else:
        dynamic_time.plot(n_s, v_s, label="Log base: " + str(b))
dynamic_time.legend(bbox_to_anchor=(0.5, 1.2), loc='upper center', borderaxespad=0., ncol=3, fancybox=True, shadow=True)
fig11.savefig("./plots/dynamic_times.png", bbox_inches='tight', pad_inches=0.2)

fig12 = plt.figure(12)
big_relative = fig12.add_subplot(111)
big_relative.grid(True)
big_relative.set_xlabel("Input Size")
big_relative.set_ylabel("Relative performance vs Nussinov\n(Higher is better)") 
for q in static_dict: 
    n_s = []
    v_s = []
    if q in [4,5,6]: 
        for n in sorted(static_dict[q]):
            n_s.append(n)
            v_s.append(static_dict[0][n]/static_dict[q][n])
        big_relative.plot(n_s, v_s, label="q: " + str(q))
for b in dynamic_dict: 
    n_s = []
    v_s = []
    if b in [4,5,6]: 
        for n in sorted(dynamic_dict[b]):
            n_s.append(n)
            v_s.append(dynamic_dict[0][n]/dynamic_dict[b][n])
        big_relative.plot(n_s, v_s, "-.", label="b: " + str(b), lw=3)
big_relative.annotate('q:6', xy=(3000, 1.93), xytext=(2900, 1.8), arrowprops=dict(facecolor='black', shrink=0.05, width=0.5, headwidth=5))
big_relative.annotate('q:5', xy=(4000, 2.39), xytext=(3900, 2.27), arrowprops=dict(facecolor='black', shrink=0.05, width=0.5, headwidth=5))
big_relative.legend(bbox_to_anchor=(0.5, 1.15), loc='upper center', borderaxespad=0., ncol=3, fancybox=True, shadow=True)
fig12.savefig("./plots/big_relative.png", bbox_inches='tight', pad_inches=0.2)