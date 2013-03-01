#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import cPickle

from datetime import datetime

data = []


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print("Usage: postproc_assimilation.py frames_dir")
        sys.exit(1)

    path = sys.argv[1]
    lst = glob.glob(os.path.join(path, "frame*"))
    N = len(lst)

    print("Loading the output data.")

    # if a pickle object exists
    pickle_path = os.path.join(path, "_pickled_cache")
    if os.path.exists(pickle_path):
        with open(pickle_path, "r") as f:
            data = cPickle.load(f)
    else:
        # open files in sequence and read the internals
        for i in range(1, N+1):
            with open(os.path.join(path, "frame%d" % i), "r") as f:
                data.append(eval(f.read()))

        # pickle to cache future use
        with open(pickle_path, "w") as f:
            cPickle.dump(data, f)

    print("Constructing time plots.")

    # extract the betas and model times
    beta = np.zeros((N, np.prod(data[0]['kriging_beta'].shape)))
    maes = np.zeros((N, 2))
    mt = []
    for i in range(N):
        mt.append(data[i]['mt'])
        beta[i,:] = data[i]['kriging_beta'][0,:]
        maes[i,0] = data[i]['model_raws_mae']
        maes[i,1] = data[i]['model_na_raws_mae']

    # plot the betas
    plt.figure()
    Nc  = beta.shape[1]
    for i in range(Nc):
        plt.subplot(Nc, 1, i+1)
        plt.plot(mt, beta[:,i])
        plt.ylabel('$\\beta_%d $ [-]' % (i+1))

    plt.xlabel('Time [-]')
    plt.savefig(os.path.join(path, "kriging_betas.png"))

    # plot maes
    plt.figure()
    plt.plot(mt, maes)
    plt.legend(['assim', 'No assim'])
    plt.ylabel('Mean abs difference [g]')
    plt.xlabel('Time [-]')
    plt.savefig(os.path.join(path, "model_maes.png"))

    print("Generating observation match for %d time points." % N)

    # plot the raws, model values at obs points and kriging
    plt.figure()
    for i in range(N):
        plt.clf()
        di = data[i]
        obs = di['kriging_obs']
        ngp = di['kriging_obs_ngp']
        m = [di['fm10_model_state'][p] for p in ngp]
        m_na = [di['fm10_model_na_state'][p] for p in ngp]
        kf = [di['kriging_field'][p] for p in ngp]
        plt.plot(m, 'ro')
        plt.plot(m_na, 'go')
        plt.plot(obs, 'bo')
        plt.plot(kf, 'ko')
        plt.legend(['assim', 'no assim', 'raws', 'kriging'])
        plt.ylabel('Fuel moisture [g]')
        plt.savefig(os.path.join(path, "image_%03d.png" % i))


    # gather all station ids
    sids = {}
    for i in range(N):
        sobs = data[i]['kriging_obs_station_ids']
        sngp = data[i]['kriging_obs_ngp']
        for sid,ngp in zip(sobs, sngp):
            if sid not in sids:
                sids[sid] = ngp

    print("Generating station plots for following stations:")
    print(sids.keys())

    plt.figure()
    for sid,ngp in sids.iteritems():
        plt.clf()

        # gather all observations from this station
        obs = []
        for j in range(N):
            dj = data[j]['kriging_obs_station_ids']
            if sid in dj:
                ndx = dj.index(sid)
                obs.append(data[j]['kriging_obs'][ndx])
            else:
                obs.append(float("nan"))

        m = [ data[j]['fm10_model_state'][ngp] for j in range(N) ]
        mna = [ data[j]['fm10_model_na_state'][ngp] for j in range(N) ]
        kfo = [ data[j]['kriging_field'][ngp] for j in range(N) ]
        plt.plot(mt, m, 'r-')
        plt.plot(mt, mna, 'g-')
        plt.plot(mt, obs, 'bo')
        plt.plot(mt, kfo, 'ko')
        plt.legend(['assim', 'no assim', 'raws', 'kriging'])
        plt.ylabel('Fuel moisture [g]')
        plt.savefig(os.path.join(path, "station_%s.png" % sid))
