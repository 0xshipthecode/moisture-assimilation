#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import cPickle
from multiprocessing import Pool, Queue, Process

from datetime import datetime

num_workers = 16
data = []


def file_loader_slave(x):
    i, fname = x

    # read in the file
    with open(fname, "r") as f:
        d = eval(f.read())

    return i, d


def plot_maker(jobs):
    # reimport what we need for plotting
    import matplotlib.pyplot as plt

    plt.figure(figsize = (12,8))

    while True:
        # retrieve next assignment (or None if end of queue)
        tmp = jobs.get()
        if tmp == None:
            break

        # decode
        (m, m_a, m_na, obs, kf, sids_pos, sids_list, fm_max, fname) = tmp
        
        plt.clf()
        plt.plot(m, 'ro')
        plt.plot(m_a, 'rs')
        plt.plot(m_na, 'go')
        plt.plot(obs, 'bo')
        plt.plot(kf, 'ko')
        plt.xticks(np.arange(len(sids_pos)), sids_list, rotation = 90)
        plt.legend(['pre-assim', 'post-assim', 'no assim', 'raws', 'kriging'])
        plt.ylabel('Fuel moisture [-]')
        plt.xlim(-0.5, len(sids_pos) - 0.5)
        plt.ylim(0.0, fm_max)
        plt.savefig(fname)


def field_plot_maker(jobs):

    #import what we need for plotting
    import matplotlib.pyplot as plt
    plt.figure(figsize = (12,8))
    
    while True:
        # retrieve next assignment (or None if end of queue)
        tmp = jobs.get()
        if tmp == None:
            break

        # decode
        (f, title, fm_max, fname) = tmp
        
        plt.clf()
        plt.imshow(f)
        plt.clim(0.0, fm_max)
        plt.title(title, fontsize = 16)
        plt.colorbar()
        plt.savefig(fname)



def plot_variogram(dists, sqdiffs, bins, title, savename):
    """
    Plots the various (halved) squared differences as a point plot and overlays
    a curve, which is the empirical variogram estimate.
    """
    bin_width = bins[1] - bins[0]

    # compute an empirical variogram
    bsqdiffs = []
    for lim in bins:
        ndx = [i for i in range(len(dists)) if (lim - bin_width <= dists[i]) and (dists[i] < lim)]
        bsqdiffs.append(np.mean([sqdiffs[i] for i in ndx]))

    # construct a plot at this time
    plt.clf()
    plt.plot(dists, sqdiffs, 'ro', markersize = 5)
    plt.plot(np.arange(bin_width, max_dist+bin_width, bin_width) - bin_width/2.0, bsqdiffs, 'b-', linewidth = 2.0)
    plt.xlabel('Distance [km]')
    plt.ylabel('Squared obs. diff [-]')
    plt.title('Variogram at time %s' % str(t_now))
    plt.savefig(savename)




def make_subdirs(path):
    """
    Construct subdirectories to store various parts of the results.
    """
    for subp in ["fm10_assim", "fm10_na", "stats", "station_xsection", "stations", "res_variograms"]:
        try:
            os.mkdir(os.path.join(path, subp))
        except OSError:
            pass


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print("Usage: postproc_assimilation.py frames_dir")
        sys.exit(1)

    path = sys.argv[1]
    what = sys.argv[2:]
    lst = glob.glob(os.path.join(path, "frame*"))
    N = len(lst)

    # construct all subdirectories
    make_subdirs(path)

    # load the data in parallel
    print("Loading the output data.")
    pool = Pool(16)
    data = [None] * N
    read_job_list = [ (i, os.path.join(path, "frame%d" % i)) for i in range(1, N+1) ]
    print len(read_job_list)
    read_res = pool.map(file_loader_slave, read_job_list)
    for i,di in read_res:
        data[i-1] = di

    # make read_res available for loading
    read_res = None

    pool.close()

    # gather all station ids that provide observations
    sids = {}
    for i in range(N):
        sobs = data[i]['kriging_obs_station_ids']
        sngp = data[i]['kriging_obs_ngp']
        for sid,ngp in zip(sobs, sngp):
            if sid not in sids:
                sids[sid] = (ngp[0] - 1, ngp[1] - 1)

    # sort keys to obtain a plotting order
    sids_list = sorted(sids.keys())
    sids_ngp = [ sids[s] for s in sids_list ]
    sids_pos = {}
    for i in range(len(sids_list)):
        sids_pos[sids_list[i]] = i

    # extract the betas, model times and kriging errors
    beta = np.zeros((N, np.prod(data[0]['kriging_beta'].shape)))
    maes = np.zeros((N, 3))
    mt = []
    ks2 = np.zeros(N)
    kerrors = []
    for i in range(N):
        mt.append(data[i]['mt'])
        beta[i,:] = data[i]['kriging_beta'][0,:]
        maes[i,0] = data[i]['model_raws_mae']
        maes[i,1] = data[i]['model_raws_mae_assim']
        maes[i,2] = data[i]['model_na_raws_mae']
        ks2[i] = data[i]['kriging_sigma2_eta']

        derrs = di['kriging_errors'][0,:]
        errs_i = np.zeros(len(sids_list))
        errs_i[:] = float("nan")
        for (s, j) in zip(di['kriging_obs_station_ids'],range(Nobs)):
            errs_i[sids_pos[s]] = derrs[j]
        kerrors.append(errs_i)

    # prepare dates for x axis of time plots
    date_ndx = np.arange(0, N, N/20)
    dates = [mt[i].strftime("%m-%d %H:%M") for i in date_ndx]
    
    if "stats" in what or "all" in what:
        print("Constructing time plots.")
    
        # plot the betas
        plt.figure(figsize=(12,8))
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
        Nc = beta.shape[1]
        for i in range(Nc):
            plt.clf()
            plt.plot(beta[:,i])
            plt.ylabel('$\\beta_%d $ [-]' % (i+1))
            plt.xlabel('Time [-]')
            plt.xticks(date_ndx, dates, rotation = 90, size = 'small')
            y = plt.ylim()
            plt.ylim(0.0, y[1])
            plt.savefig(os.path.join(path, "stats", "kriging_beta_%d.png" % (i+1)))

        # plot maes
        plt.clf()
        plt.plot(maes)
        plt.legend(['forecast', 'analysis', 'no assim'])
        plt.ylabel('Mean abs difference [g]')
        plt.xlabel('Time [-]')
        plt.xticks(date_ndx, dates, rotation = 90, size = 'small')
        y = plt.ylim()
        plt.ylim(0, y[1])
        plt.savefig(os.path.join(path, "stats", "model_maes.png"))

        # plot the etas
        plt.clf()
        plt.plot(ks2)
        plt.ylabel('Microscale variability variance [-]')
#        plt.xlabel('Time [-]')
        plt.xticks(date_ndx, dates, rotation = 90, size = 'small')
        y = plt.ylim()
        plt.ylim(0, y[1])
        plt.savefig(os.path.join(path, "stats", "eta_variance.png"))
    

    # find the maximal values of the plotted variables
    fm_max = 0.0
    fm_max_na = 0.0
    fm_max_assim = 0.0
    for i in range(N):
        di = data[i]
        fm_max_assim = max(fm_max_assim, np.amax(di['fm10_model_state_assim']))
        fm_max_na = max(fm_max_na, np.amax(di['fm10_model_na_state']))
        fm_max = max(fm_max, np.amax(di['fm10_model_state']), fm_max_na, fm_max_assim,
                     np.amax(di['kriging_field']), np.amax(di['kriging_obs']))

    # artificially chop the max so that in case there are extreme values at least something is visible
    fm_max = min(fm_max, 0.6)
    fm_max_assim = min(fm_max_assim, 0.6)
    fm_max_na = min(fm_max_na, 0.6)

    plot_queue = Queue()
    field_queue = Queue()

    # plot the raws, model values at obs points and kriging
    if "fields" in what or "all" in what:
        print("Generating observation/model matchups for %d time points." % N)

        for i in range(N):
            di = data[i]
            dobs = di['kriging_obs']
            Nobs = len(dobs)

            # fill out observations we have for this frame
            obs = np.zeros(len(sids_list))
            obs[:] = float("nan")
            for (s, j) in zip(di['kriging_obs_station_ids'],range(Nobs)) :
                obs[sids_pos[s]] = dobs[j]

            m = [di['fm10_model_state'][p] for p in sids_ngp]
            m_a = [di['fm10_model_state_assim'][p] for p in sids_ngp]
            m_na = [di['fm10_model_na_state'][p] for p in sids_ngp]
            kf = [di['kriging_field'][p] for p in sids_ngp]
            fname = os.path.join(path, "image_%03d.png" % i)

            plot_queue.put((m, m_a, m_na, obs, kf, sids_pos, sids_list, fm_max, fname))

            field_queue.put((di['fm10_model_state_assim'],
                             'Fuel moisture state %d at %s (ASSIM)' % (i, str(mt[i])), fm_max_assim,
                             os.path.join(path, "fm10_assim", "fm10_assim_field_%03d.png" % i)))

            field_queue.put((di['fm10_model_na_state'],
                             'Fuel moisture state %d at %s (NA)' % (i, str(mt[i])), fm_max_na,
                             os.path.join(path, "fm10_na", "fm10_na_field_%03d.png" % i)))


        # start up the workers and process the queue
        workers = []
        for i in range(num_workers): 
            # end-of-queue marker (one for each worker)
            plot_queue.put(None)

            # create a new worker and add it to the pool
            tmp = Process(target=plot_maker, args=(plot_queue,))
            tmp.start()
            workers.append(tmp)

        # wait for workers to complete
        for worker in workers:
            worker.join()

        workers = []
        for i in range(num_workers): 
            # end-of-queue marker (one for each worker)
            field_queue.put(None)

            # create a new worker and add it to the pool
            tmp = Process(target=field_plot_maker, args=(field_queue,))
            tmp.start()
            workers.append(tmp)

        # wait for workers to complete
        for worker in workers:
            worker.join()


    # print kriging statistics if statistics requested
    if "stats" in what or "all" in what:
        plt.figure(figsize=(12,8))
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
        err_variance = []
        kerrors = []
        for i in range(N):
            err_i = data[i]['kriging_errors']
            err_variance.append(1.0 / len(err_i) * np.sum(err_i ** 2))

            derrs = di['kriging_errors'][0,:]
            errs_i = np.zeros(len(sids_list))
            errs_i[:] = float("nan")
            for (s, j) in zip(di['kriging_obs_station_ids'],range(Nobs)):
                errs_i[sids_pos[s]] = derrs[j]
            kerrors.append(errs_i)

        
        plt.plot(err_variance)
        plt.title('Variance of model error vs. time')
        plt.ylabel('Kriging Error variance [-]')
        plt.xlabel('Time')
        date_ndx = np.arange(0, N, N/20)
        dates = [mt[i].strftime("%m-%d %H:%M") for i in date_ndx]
        plt.xticks(date_ndx, dates, rotation = 90, size = 'small')
        plt.savefig(os.path.join(path, "error_variance.png"))

        plt.clf()
        kecc = np.ma.corrcoef(np.ma.array(kerrors, mask = np.isnan(kerrors)), rowvar = False)
        plt.imshow(kecc, interpolation='nearest')
        plt.colorbar()
        plt.title('Correlation coefficients of kriging errors')
        plt.clim([-1, 1])
        ax = plt.gca()
        ax.set_xticks(np.arange(len(sids_list)))
        ax.set_xticklabels(sids_list, rotation = 90)
        ax.set_yticks(np.arange(len(sids_list)))
        ax.set_yticklabels(sids_list)
        plt.savefig(os.path.join(path, "kriging_error_correlations.png"))

        plt.clf()
        trindx = np.triu_indices(kecc.shape[0], 1)
        plt.hist(kecc[trindx], int(trindx[0].shape[0]**0.5), normed = True)
        plt.xlim([-1, 1])
        plt.title("Histogram of correlation coefficients of kriging errors")
        plt.savefig(os.path.join(path, "kriging_errors_histogram.png"))


    # plot agreement with simulated and station data
    if "station_plots" in what or "all" in what:
        plt.figure(figsize = (12,8))
        for sid,ngp in sids.iteritems():
            plt.clf()
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)

            # gather all observations from this station
            obs = []
            for j in range(N):
                dj = data[j]['kriging_obs_station_ids']
                if sid in dj:
                    ndx = dj.index(sid)
                    obs.append(data[j]['kriging_obs'][ndx])
                else:
                    obs.append(float("nan"))

            m    = [ data[j]['fm10_model_state'][ngp] for j in range(N) ]
            m_a  = [ data[j]['fm10_model_state_assim'][ngp] for j in range(N) ]
            m_na = [ data[j]['fm10_model_na_state'][ngp] for j in range(N) ]
            kfo  = [ data[j]['kriging_field'][ngp] for j in range(N) ]
            kv   = [ data[j]['kriging_variance'][ngp] for j in range(N) ]
            mv   = [ data[j]['fm10_model_var'][ngp] for j in range(N) ]
            plt.plot(m, 'r--')
            plt.plot(m_a, 'r-')
            plt.plot(m_na, 'g-')
            plt.plot(obs, 'bx-')
            plt.plot(kfo, 'k+')
            date_ndx = np.arange(0, N, N/20)
            dates = [mt[i].strftime("%m-%d %H:%M") for i in date_ndx]
            plt.xticks(date_ndx, dates, rotation = 90, size = 'small')
            plt.legend(['pre-assim', 'post-assim', 'no assim', 'raws', 'kriging'])
            plt.ylabel('Fuel moisture [g]')
            plt.savefig(os.path.join(path, "station_%s.png" % sid))

            plt.clf()
            plt.plot(np.log10(kv), 'ko')
            plt.plot(np.log10(mv), 'ro')
            plt.xticks(date_ndx, dates, rotation = 90, size = 'small')
            plt.legend([ 'kriging var', 'model var' ])
            plt.ylabel('Log10 variance')
            plt.savefig(os.path.join(path, "station_%s_var.png" % sid))



    if "variograms" in what or "all" in what:
        plt.figure(figsize = (12,8))

        dists = []
        sqdiffs = []

        # we must load stations pertinent to the simulation (problem)
        # load station data from files
        with open(cfg['station_list_file'], 'r') as f:
            si_list = f.read().split('\n')
        si_list = filter(lambda x: len(x) > 0 and x[0] != '#', map(string.strip, si_list))

        stations = {}
        for code in si_list:
            mws = MesoWestStation(code)
            mws.load_station_info(os.path.join(cfg["station_info_dir"], "%s.info" % code))
            stations[code] = mws

        # first construct a distance matrix for all stations
        st_dists = np.zeros((len(sids_list), len(sids_list)))
        for j in range(len(sids_list)):
            lonj, latj = j.get_position()
            for k in range(j+1, len(sids_list)):
                lonk, latk = k.get_position()
                 if j in sids_pos and k in sids_pos:
                     d = great_circle_distance(lonj, latj, lonk, latk)
                     st_dists[sids_pos[j], sids_pos[k]] = d
                     st_dists[sids_pos[k], sids_pos[j]] = d
                             

        # accumulate data over all time points
        for i in range(N):

            # retrieve valid measurements, rest is nan
            di = data[i]
            Nobs = len(di['kriging_obs'])
            obs = np.zeros(len(sids_list))
            obs[:] = float("nan")
            for (s, j) in zip(di['kriging_obs_station_ids'],range(Nobs)) :
                obs[sids_pos[s]] = dobs[j]

            # loop through pairs
            for j in range(Nobs):
                for k in range(j+1, Nobs):
                    if not np.isnan(obs[j] + obs[k]):
                        dists.append(st_dists[j,k])
                        sq_diffs.append(0.5 * (obs[j] - obs[k])**2)


        plot_variogram(dists, sqdiffs, np.arange(20, 500, 20), 'Variogram from TSM residuals',
                       os.path.join(path, "res_variograms", "kriging_residual_variogram.png"))
