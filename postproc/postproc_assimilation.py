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


def file_loader(x):
    i, fname = x
    with open(fname, "r") as f:
        d = eval(f.read())
    return (i, d)


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


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print("Usage: postproc_assimilation.py frames_dir")
        sys.exit(1)

    path = sys.argv[1]
    lst = glob.glob(os.path.join(path, "frame*"))
    N = len(lst)

    print("Loading the output data.")

    pool = Pool(16)

    # if a pickle object exists
#    pickle_path = os.path.join(path, "_pickled_cache")
#    if os.path.exists(pickle_path):
#        with open(pickle_path, "r") as f:
#            data = cPickle.load(f)
#    else:
        # open files in sequence and read the internals
    # create a set of jobs to read in the files
    data = [None] * N
    read_job = [ (i, os.path.join(path, "frame%d" % i)) for i in range(1, N+1) ]
    read_res = pool.map(file_loader, read_job)
    for i,di in read_res:
        data[i-1] = di

    pool.close()

#    for i in range(1, N+1):
#        with open(os.path.join(path, "frame%d" % i), "r") as f:
#            data.append(eval(f.read()))

        # pickle to cache future use
#        with open(pickle_path, "w") as f:
#            cPickle.dump(data, f)

    print("Constructing time plots.")

    # extract the betas and model times
    beta = np.zeros((N, np.prod(data[0]['kriging_beta'].shape)))
    maes = np.zeros((N, 3))
    mt = []
    ks2 = np.zeros(N)
    for i in range(N):
        mt.append(data[i]['mt'])
        beta[i,:] = data[i]['kriging_beta'][0,:]
        maes[i,0] = data[i]['model_raws_mae']
        maes[i,1] = data[i]['model_raws_mae_assim']
        maes[i,2] = data[i]['model_na_raws_mae']
        ks2[i] = data[i]['kriging_sigma2_eta']

    # prepare dates for x axis of time plots
    date_ndx = np.arange(0, N, N/20)
    dates = [mt[i].strftime("%m-%d %H:%M") for i in date_ndx]

    # plot the betas
    plt.figure(figsize=(12,8))
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
    Nc  = beta.shape[1]
    for i in range(Nc):
        plt.clf()
        plt.plot(beta[:,i])
        plt.ylabel('$\\beta_%d $ [-]' % (i+1))
        plt.xlabel('Time [-]')
        plt.xticks(date_ndx, dates, rotation = 90, size = 'small')
        plt.savefig(os.path.join(path, "kriging_beta_%d.png" % (i+1)))

    # plot maes
    plt.clf()
    plt.plot(maes)
    plt.legend(['forecast', 'analysis', 'no assim'])
    plt.ylabel('Mean abs difference [g]')
    plt.xlabel('Time [-]')
    plt.xticks(date_ndx, dates, rotation = 90, size = 'small')
    plt.savefig(os.path.join(path, "model_maes.png"))

    # plot the etas
    plt.clf()
    plt.plot(ks2)
    plt.ylabel('Kriging eta variance [-]')
    plt.xlabel('Time [-]')
    plt.xticks(date_ndx, dates, rotation = 90, size = 'small')
    plt.savefig(os.path.join(path, "eta_variance.png"))
    
    

    print("Generating observation/model matchups for %d time points." % N)

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

    print("Generating station/image plots for following stations:")
    print(sids_list)

    # find the maximal values of the plotted variables
    fm_max = 0.0
    for i in range(N):
        di = data[i]
        fm_max = max(fm_max, np.amax(di['fm10_model_state']), np.amax(di['fm10_model_na_state']),
                     np.amax(di['kriging_field']), np.amax(di['kriging_obs']))

    # artificially chop the max so that in case there are extreme values at least something is visible
    fm_max = min(fm_max, 0.6)

    plot_queue = Queue()

    # plot the raws, model values at obs points and kriging
#    plt.figure(figsize = (12,8))
    for i in range(N):
#        plt.clf()
        di = data[i]
        dobs = di['kriging_obs']
        Nobs = len(dobs)

        # fill out observations we have for this frame
        obs = np.zeros(len(sids_list))
        obs[:] = float("nan")
        for (s,j) in zip(di['kriging_obs_station_ids'],range(Nobs)) :
            obs[sids_pos[s]] = dobs[j]

        m = [di['fm10_model_state'][p] for p in sids_ngp]
        m_a = [di['fm10_model_state_assim'][p] for p in sids_ngp]
        m_na = [di['fm10_model_na_state'][p] for p in sids_ngp]
        kf = [di['kriging_field'][p] for p in sids_ngp]
        fname = os.path.join(path, "image_%03d.png" % i)

        plot_queue.put((m, m_a, m_na, obs, kf, sids_pos, sids_list, fm_max, fname))

        # plt.plot(m, 'ro')
        # plt.plot(m_a, 'rs')
        # plt.plot(m_na, 'go')
        # plt.plot(obs, 'bo')
        # plt.plot(kf, 'ko')
        # plt.xticks(np.arange(len(sids_pos)), sids_list, rotation = 90)
        # plt.legend(['pre-assim', 'post-assim', 'no assim', 'raws', 'kriging'])
        # plt.ylabel('Fuel moisture [g]')
        # plt.ylim(0.0, fm_max)
        # plt.savefig(os.path.join(path, "image_%03d.png" % i))


    plt.figure(figsize=(12,8))
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
    err_variance = []
    for i in range(N):
        err_i = data[i]['kriging_errors']
        err_variance.append(1.0 / len(err_i) * np.sum(err_i ** 2))
        
    plt.plot(err_variance)
    plt.title('Variance of model error vs. time')
    plt.ylabel('Kriging Error variance [-]')
    plt.xlabel('Time')
    date_ndx = np.arange(0, N, N/20)
    dates = [mt[i].strftime("%m-%d %H:%M") for i in date_ndx]
    plt.xticks(date_ndx, dates, rotation = 90, size = 'small')
    plt.savefig(os.path.join(path, "error_variance.png"))
        

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

