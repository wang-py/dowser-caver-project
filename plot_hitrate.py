from dowser_hitrate_functions import *

if __name__ == '__main__':
    exp_waters = sys.argv[1]
    dowser_waters = sys.argv[2:]
    exp_water_arr = read_exp_water(exp_waters)
    hitrate = []
    energy_cutoff = np.linspace(0, -12, len(dowser_waters))
    for one_dowser in dowser_waters:
        dowser_water_arr = dowser_water_pruning(one_dowser)
        hitrate.append(hit_detection(exp_water_arr, dowser_water_arr, 1.4))
    hitrate = np.array(hitrate)
    plot_hitrate_vs_cutoff(hitrate, energy_cutoff)
    pass