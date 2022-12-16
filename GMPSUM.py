import numpy as np
import copy
import math
import timeit
import sys

def calculationHProb(m,sigSize):
    P = np.zeros((m,m))
    for q in range(m):
        for k in range(m):
            if k == 0:
                P[k][q] = 1
            else:
                if k > q:
                    P[k][q] = 0
                else:
                    P[k][q] = ((1/sigSize)*P[k-1][q-1]) + (((sigSize - 1)/sigSize) * P[k][q-1])
    return P
def cleanString(s):
    for i in range(10):
        s = s.replace(str(i), '')
    lenS = len(s)
    for i in range(lenS):
        s = s.strip()
    return s
def getStrings(filename, n, lenalphabet, withspace = False):
    if lenalphabet == 100:
        f = open(filename, 'r', encoding='iso-8859-15')
    else:
        f = open(filename, 'r')
    f.readline()
    _strings = []
    for i in range(n):
        _strings.append(cleanString(f.readline()))
        if withspace == True:
            f.readline()
    return _strings
def init(alphabetSize = 4):
    beta = 50
    k_best = 7
    if alphabetSize == 4:
        alphabet = ['A', 'C', 'G', 'T']
    elif alphabetSize == 20:
        alphabet = ['A', 'C', 'D','E','F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'S', 'T', 'V', 'W', 'X', 'Y']
    elif alphabetSize == 2:
        alphabet = ['a', 'b']
    elif alphabetSize == 8:
        alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
    elif alphabetSize == 24:
        alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j','k', 'l', 'm', 'n', 'o','p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x']
    else:
        print("Error. alphabet is not defined for this dataset!")
    alpha = 30 #this parameter should be defined using a function n = 60 ==> alpha = 40
    return alphabet, beta, k_best, alpha
def feasible(locations, alphabet, S, lenstrings):
    for i in range(lenstrings):
        if S[i][locations[i]:].find(alphabet) == -1:
            return False
    return True
def successor(locations, alphabet, S, lenstrings):
    locs = []
    for i in range(lenstrings):
        locs.append(locations[i] + S[i][locations[i]:].find(alphabet) + 1)
    return locs
def calculateUpperBound(S, alphabet, lens):
    _ub = 0
    lenS = len(S)
    lenalphabet = len(alphabet)
    for j in range(lenalphabet):
        temp_ub = 10000
        for i in range(lenS):
            if S[i][len(S[i]) - lens[i]:].count(alphabet[j]) < temp_ub:
                temp_ub = S[i][len(S[i]) - lens[i]:].count(alphabet[j])
        _ub = _ub + temp_ub
    return _ub
def muCalculation(x):
    p = 1
    for i in range(len(x)):
        if len(x) > 100:
            p = p * (x[i] / 4)
        if len(x) <= 100:
            p = p * x[i]
    p = p ** (1/len(x))
    return p
def stdCalculation(x,mu_):
    sum_  = 0
    for i in range(len(x)):
        if len(x) <= 100:
            t = np.log(x[i]/mu_)
            t = t ** 2
            sum_ += t
        if (x[i] / mu_) != 0 and len(x) > 100:
            t = np.log(x[i]/mu_)
            t = t ** 2
            sum_ += t
    sum_ = sum_ / len(x)
    sum_ = math.sqrt(sum_)
    sigma = math.e ** sum_
    return sigma
def ocurranceOfStrings(S, lens, a):
    x = []
    for i in range(len(S)):
        x.append(S[i][len(S[i]) - lens[i]:].count(a))
    return x
def calculateScore(C, S, lambda_, alphabet, lenstrings, P):
    lens = []
    for i in range(lenstrings):
        lens.append(len(S[i]) - C[i])
    gm = 0.0
    ub = calculateUpperBound(S, alphabet, lens)
    if ub == 0:
        return 0
    for i in range(len(alphabet)):
        ca = ocurranceOfStrings(S, lens, alphabet[i])
        mu_ = muCalculation(ca)
        if mu_ == 0:
            continue
        std_ = stdCalculation(ca, mu_)
        h_temp = (mu_ / std_) * (min(ca) / ub)
        gm = gm + h_temp
    l_max = min(lens)
    psum = 0.0
    for k in range(l_max):
        h_temp = 1.0
        for i in range(lenstrings):
            h_temp = h_temp * P[k][lens[i]]
        psum = psum + h_temp
    h_value =  ((lambda_ * gm) + (1.0 - lambda_) * psum)
    return h_value
def checkDominance(c, best):
    lenc = len(c)
    for i in range(lenc):
        if c[i] <= best[i]:
            return False
    return True
def kBestListFilter(C,h,kbest):
    kbestlist = []
    htemp = copy.deepcopy(h)
    for i in range(kbest):
        t = -1
        ind = 0
        for j in range(len(h)):
            if h[j] > t:
                ind = j
                t = h[j]
        h[ind] = -1
        kbestlist.append(C[ind])
    lenC = len(C)
    Cnewlist = []
    hnewlist = []
    counter = 0
    for i in range(lenC):
        flag = True
        for j in range(kbest):
            if checkDominance(C[i], kbestlist[j]) == True:
                flag = False
                counter += 1
                break
        if flag == True:
            Cnewlist.append(C[i])
            hnewlist.append(htemp[i])
    return Cnewlist, hnewlist, counter
def keepBetaBest(C,h,Beta):
    _rtn = []
    for i in range(Beta):
        t = -1
        ind = 0
        for j in range(len(h)):
            if h[j] > t:
                ind = j
                t = h[j]
        h[ind] = -1
        _rtn.append(C[ind])
    return _rtn
def findKMin(C,S):
    lens = []
    lenstrings = len(S)
    for i in range(lenstrings):
        lens.append(len(S[i]) - C[i])
    _rtn = min(lens)
    return _rtn
def main(lns, lna, dta, inalphabet, flg, Beta, lmbd = 0.0):
    lenalphabet = lna
    lenstrings = lns
    dataset = dta
    filename = dataset
    alphabet = inalphabet
    lambda_ = lmbd
    S = getStrings(filename, lenstrings, lenalphabet, flg) #ST
    beta = Beta
    tls = []
    for i in range(len(S)):
        tls.append(len(S[i]))
    tls = max(tls)
    P = calculationHProb(tls, lenalphabet)
    B = []
    B.append(np.zeros(lenstrings, dtype= 'int'))
    lenlcs = 0
    starttime = timeit.default_timer()
    les = []
    for i in range(lenstrings):
        les.append(len(S[i]))
    while True:
        lenB = len(B)
        C = []
        h = []
        for i in range(lenB):
            for j in range(lenalphabet):
                if feasible(B[i], alphabet[j], S, lenstrings) == True:
                    C.append(successor(B[i], alphabet[j], S, lenstrings))
        lenC = len(C)
        for i in range(lenC):
            h.append(calculateScore(C[i], S, lambda_, alphabet, lenstrings, P))
        if lenC == 0:
            break
        else:
            #C, h, co = kBestListFilter(C, h, k_best)
            if lenC > beta:
                B = keepBetaBest(C,h, beta)
            else:
                B = C
        lenlcs +=1
    stoptime = timeit.default_timer()
    return (lenlcs), stoptime - starttime

def GMPSUM(name, _Beta, lmbd):
    if name == "ACO-Random":
        lenalphabet = [4, 20]
        lenstrings = [10, 15, 20, 25, 40,60,80,100,150,200]
        datasets = ["rnd"]
        for i in range(len(datasets)):
            for j in range(len(lenalphabet)):
                for k in range(len(lenstrings)):
                    if lenalphabet[j] == 4:
                        inalphabet = ['A', 'C', 'G', 'T']
                    elif lenalphabet[j] == 20:
                        inalphabet = ['A', 'C', 'D','E','F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'S', 'T', 'V', 'W', 'X', 'Y']
                    fname = str(len(inalphabet)) + "_" + str(lenstrings[k]) + "_600." + datasets[i]
                    lenlcs, tim = main(lenstrings[k], lenalphabet[j], fname, inalphabet, False, _Beta, lmbd)
                    print(fname, lenlcs, tim)
    if name == "ACO-Rat":
        lenalphabet = [4, 20]
        lenstrings = [10, 15, 20, 25, 40,60,80,100,150,200]
        datasets = ["rat"]
        for i in range(len(datasets)):
            for j in range(len(lenalphabet)):
                for k in range(len(lenstrings)):
                    if lenalphabet[j] == 4:
                        inalphabet = ['A', 'C', 'G', 'T']
                    elif lenalphabet[j] == 20:
                        inalphabet = ['A', 'C', 'D','E','F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'S', 'T', 'V', 'W', 'X', 'Y']
                    fname = str(len(inalphabet)) + "_" + str(lenstrings[k]) + "_600." + datasets[i]
                    lenlcs, tim = main(lenstrings[k], lenalphabet[j], fname, inalphabet, False, _Beta, lmbd)
                    print(fname, lenlcs, tim)
    if name == "ACO-Virus":
        lenalphabet = [4, 20]
        lenstrings = [10, 15, 20, 25, 40,60,80,100,150,200]
        datasets = ["virus"]
        for i in range(len(datasets)):
            for j in range(len(lenalphabet)):
                for k in range(len(lenstrings)):
                    if lenalphabet[j] == 4:
                        inalphabet = ['A', 'C', 'G', 'T']
                    elif lenalphabet[j] == 20:
                        inalphabet = ['A', 'C', 'D','E','F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'S', 'T', 'V', 'W', 'X', 'Y']
                    fname= str(len(inalphabet)) + "_" + str(lenstrings[k]) + "_600." + datasets[i]
                    lenlcs, tim  = main(lenstrings[k], lenalphabet[j], fname, inalphabet, False, _Beta, lmbd)
                    print(fname, lenlcs, tim)
    if name == "BB":
        lenalphabet = [2, 4, 8, 24]
        lenstrings = [10, 100]
        for i in range(len(lenalphabet)):
            if lenalphabet[i] == 2:
                inalphabet = ['a', 'b']
            if lenalphabet[i] == 4:
                inalphabet = ['a', 'b', 'c', 'd']
            if lenalphabet[i] == 8:
                inalphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
            if lenalphabet[i] == 24:
                inalphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j','k', 'l', 'm', 'n', 'o','p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x']
            for j in range(len(lenstrings)):
                answers = []
                _timav = []
                for k in range(10):
                    fname = str(lenalphabet[i]) + "_" + str(lenstrings[j]) + "_1000.het0.1." + str(k + 1)
                    x, tim = main(lenstrings[j], lenalphabet[i], fname, inalphabet, False, _Beta, lmbd)
                    answers.append(x)
                    _timav.append(tim)
                print(str(lenalphabet[i]) + "_" + str(lenstrings[j]) + "_1000.het0.1", np.mean(answers), np.mean(_timav))

    if name == "SARS-CoV-2":
        inalphabet = ['A', 'C', 'G', 'T']
        lenalphabet = len(inalphabet)
        lenstrings = [10,20,30,40,50,60,70,80,90,100,110]
        for i in range(len(lenstrings)):
            fname = 'cov_' + str(lenstrings[i]) + "_" + str(lenalphabet)
            lenlcs, tim = main(lenstrings[k], inalphabet, fname, inalphabet, False, _Beta, lmbd)
            print(fname, lenlcs, tim)

    if name == "ES":
        lenalphabet = [2, 10, 25, 100]
        lenstrings = [10, 50, 100]
        for i in range(len(lenalphabet)):
            if lenalphabet[i] == 2:
                inalphabet = ['0', '1']
            if lenalphabet[i] == 10:
                inalphabet = ['&', '(', '#', '%', '\"', '\'', '!', ')', '$', '*']
            if lenalphabet[i] == 25:
                inalphabet = ['"', '8', '%', '$', '0', '&', '9', ')', '1', '*', '7', '!', '5', '+', '6', '(', ',', '3', '-', '/', '4', '.', "'", '2', '#']
            if lenalphabet[i] == 100:
                inalphabet = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
                            '#', '+', ',', '*', '!', '/', '\"', '.', '(', ')',
                            '$', '&', '-', '%', '\'', 'a', 'b', 'c', 'd', 'e',
                            'f', 'g', 'h', 'i', 'g', 'k', 'l', 'm', 'n', 'o',
                            'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y',
                            'z', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I',
                            'G', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S',
                            'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '~', '`', ';',
                            '}', '^', 'þ', '€', '=', 'ƒ', '>', '@', '<', ']'
                            '|', ':', '\\', '?', '_', '[', '{', '', '', '', '']
            for j in range(len(lenstrings)):
                answers = []
                _timav = []
                for k in range(50):
                    fname = "ES_" + str(lenstrings[j]) + "_" + str(lenalphabet[i]) + "_" + str(k+1) + ".txt"
                    lenlcs, tim = main(lenstrings[j], lenalphabet[i], fname, inalphabet, False, _Beta, lmbd)
                    answers.append(lenlcs)
                    _timav.append(tim)
                print("ES_" + str(lenstrings[j]) + "_" + str(lenalphabet[i]), np.mean(answers), np.mean(_timav))

if __name__ == "__main__":
    datasetName = sys.argv[1]
    beamWidth = int(sys.argv[2])
    lmbd = float(sys.argv[3])
    GMPSUM(datasetName, beamWidth, lmbd)