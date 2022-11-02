from cmath import inf
import numpy as np
import copy
import timeit
import random
import sys
import os

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

def cleanString(s, _Num):
    if _Num == False:
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
    _Num = False
    if "ES" in filename:
        _Num = True
    for i in range(n):
        _strings.append(cleanString(f.readline(), _Num))
        if withspace == True:
            f.readline()
    return _strings

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
def k_analytic_score(C, S, _k, lenstrings, P):
    lens = []
    for i in range(lenstrings):
        lens.append(len(S[i]) - C[i])
    h_value = 1.0
    for i in range(lenstrings):
        h_value = h_value * P[_k][lens[i]]
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
def allmins(C,S):
    lens = []
    lenstrings = len(S)
    for i in range(lenstrings):
        lens.append(len(S[i]) - C[i])
    return lens

def approximate_k(sigSize, l, pi):
    term = 1.0
    if pi < 1e-250:
        return 0
    pxInterm = np.log(pi) + l * np.log(sigSize)
    if pxInterm >= 300:
        term = 1.0
    else:
        px = np.exp(pxInterm)
        tildaP = px * (-1.0 - pi/2.0)
        if abs(tildaP) < 1e-15:
            return -tildaP
        term = np.exp(tildaP)
    return 1.0 - term

def stable_eq(sigSize, l, x):
    step = 20
    val = x
    powe = l
    while powe > x:
        val = val ** (sigSize ** step)
        powe = powe - step
    val = val ** (sigSize ** powe)
    return 1.0 - val

def sf(C, S, alphabet, lenstrings, P, k):
    lens = calculate_lens(lenstrings, C, S)
    pi = 1.0
    for i in range(lenstrings):
        pi *= P[k][lens[i]]
    if pi >= 1e-10:
        pi = 1.0 - pi
        pi = stable_eq(len(alphabet),k , pi)
    else:
        pi = approximate_k(len(alphabet),k ,pi)
    return pi

def ex_score(C, S, alphabet, lenstrings, P, l_min, l_max):
    if l_min > l_max or l_min <0:
        return 0
    if l_min == l_max:
        return sf(C,S,alphabet, lenstrings,P, l_min)
    mid = int((l_min + l_max) / 2)
    mid_val = sf(C,S,alphabet, lenstrings ,P, mid)
    if mid_val <= 1e-6:
        return (l_min - mid + 1) * mid_val + ex_score(C,S,alphabet, lenstrings, P, l_min, mid - 1)
    if mid_val >= (1 - 1e-6):
        return mid - l_min + 1 + ex_score(C, S,alphabet, lenstrings, P, mid + 1, l_max)
    else:
        return mid_val + ex_score(C, S, alphabet, lenstrings, P, l_min, mid - 1) + ex_score(C, S, alphabet, lenstrings, P, mid + 1, l_max)

def calculate_lens(lenstrings, C, S):
    lens = []
    for i in range(lenstrings):
        lens.append(len(S[i]) - C[i])
    return lens

def GCoV_score(C, S, alphabet, lenstrings, gamma):
    lens = calculate_lens(lenstrings, C, S)
    ub = calculateUpperBound(S, alphabet, lens)
    Mean = np.mean(lens)
    Std = np.std(lens)
    gamma = (0.0036 * lenstrings) - 0.0161
    if Std == 0:
        Std = 1e-15
    h_value = (((ub ** (1/2))) * (Mean ** 2)) / (Std ** (gamma))
    return h_value

def ubexaminer(S, alphabet):
    lens = []
    for i in range(len(S)):
        lens.append(len(S[i]))
    ltlens = max(lens)
    ub = calculateUpperBound(S, alphabet, lens) / ltlens
    return ub

def LCSubStr(X, Y):
    m = len(X)
    n = len(Y)
    LCSuff = [[0 for k in range(n+1)] for l in range(m+1)]
    result = 0
    for i in range(m + 1):
        for j in range(n + 1):
            if (i == 0 or j == 0):
                LCSuff[i][j] = 0
            elif (X[i-1] == Y[j-1]):
                LCSuff[i][j] = LCSuff[i-1][j-1] + 1
                result = max(result, LCSuff[i][j])
            else:
                LCSuff[i][j] = 0
    return result

def S2D(S):
    lenstrings = len(S)
    lens = []
    for i in range(len(S)):
        lens.append(len(S[i]))
    ltlens = min(lens)
    n = int(lenstrings / 2)
    alpha = 20
    ans = []
    for i in range(n):
        xind = random.randint(0, ltlens - alpha)
        xnst1 = random.randint(0, lenstrings - 1)
        xnst2 = random.randint(0, lenstrings - 1)
        while xnst1 == xnst2:
            xnst2 = random.randint(0, lenstrings - 1)
        lensub = LCSubStr(S[xnst1][xind:(xind + alpha)], S[xnst2][xind:(xind + alpha)])
        ans.append(lensub)
    a = np.mean(ans)
    if a <= 5.0:
        return "Uncorrelated"
    else:
        return "Correlated"
def main(lns, lna, dta, inalphabet, flg, _Beta):
    lenalphabet = lna
    lenstrings = lns
    dataset = dta
    filename = dataset
    alphabet = inalphabet
    beta = _Beta
    S = getStrings(filename, lenstrings, lenalphabet, flg)
    _Alg = ""
    _Type = S2D(S)
    if _Type == "Correlated":
        _Alg = "k_cor"
    else:
        _Ubs = ubexaminer(S, alphabet)
        if _Ubs >= 0.9:
            _Alg = "BS-Ex"
        elif _Ubs < 0.9 and _Ubs >= 0.53:
            _Alg = "k_uncor"
        else:
            _Alg = "GCoV"
    _T = []
    for i in range(len(S)):
        _T.append(len(S[i]))
    _Lmin = max(_T)
    P = calculationHProb(_Lmin, lenalphabet)
    B = []
    B.append(np.zeros(lenstrings, dtype= 'int'))
    lenlcs = 0
    starttime = timeit.default_timer()


    if _Alg == "k_cor":
        lenlcs = 0
        while True:
            lenB = len(B)
            C = []
            h = []
            for i in range(lenB):
                for j in range(lenalphabet):
                    if feasible(B[i], alphabet[j], S, lenstrings) == True:
                        C.append(successor(B[i], alphabet[j], S, lenstrings))
            lenC = len(C)
            _K = 1000000
            for i in range(lenC):
                q_min = findKMin(C[i],S)
                if q_min < _K:
                    _K = q_min
            _K = _K - 31
            _K = int((_K) / lenalphabet)
            if _K <= 0:
                _K = 1
            for i in range(lenC):
                h.append(k_analytic_score(C[i], S, _K, lenstrings, P))
            if lenC == 0:
                break
            else:
                if lenC > beta:
                    B = keepBetaBest(C,h, beta)
                else:
                    B = C
            lenlcs +=1
        stoptime = timeit.default_timer()
        return lenlcs, (stoptime - starttime)
    
    elif _Alg == "k_uncor":
        lenlcs = 0
        while True:
            lenB = len(B)
            C = []
            h = []
            for i in range(lenB):
                for j in range(lenalphabet):
                    if feasible(B[i], alphabet[j], S, lenstrings) == True:
                        C.append(successor(B[i], alphabet[j], S, lenstrings))
            lenC = len(C)
            _K = 0
            for i in range(lenC):
                q_min = findKMin(C[i],S)
                if q_min > _K:
                    _K = q_min
            _T = 1.8233 - (0.1588 * np.log(lenstrings)) 
            _K = _K * _T
            _K = int((_K) / lenalphabet)
            if _K <= 0:
                _K = 1
            for i in range(lenC):
                h.append(k_analytic_score(C[i], S, _K, lenstrings, P))
            if lenC == 0:
                break
            else:
                if lenC > beta:
                    B = keepBetaBest(C,h, beta)
                else:
                    B = C
            lenlcs +=1
        stoptime = timeit.default_timer()
        return lenlcs, (stoptime - starttime)

    elif _Alg == "BS-Ex":
        lenlcs = 0
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
                lenss = calculate_lens(lenstrings, C[i], S)
                l_maxx = min(lenss)
                ex = ex_score(C[i], S, alphabet, lenstrings, P, 0, l_maxx)
                h.append(ex)
            if lenC == 0:
                break
            else:
                if lenC > beta:
                    B = keepBetaBest(C,h, beta)
                else:
                    B = C
            lenlcs +=1
        stoptime = timeit.default_timer()
        return lenlcs, (stoptime - starttime)
    elif _Alg == "GCoV":
        lenlcs = 0
        gamma = (0.0036 * lenstrings) - 0.0161
        while True:
            lenB = len(B)
            C = []
            h = []
            mainS = copy.deepcopy(S)
            for i in range(lenB):
                for j in range(lenalphabet):
                    if feasible(B[i], alphabet[j], S, lenstrings) == True:
                        C.append(successor(B[i], alphabet[j], S, lenstrings))
            lenC = len(C)
            for i in range(lenC):
                h.append(GCoV_score(C[i], S, alphabet, lenstrings, gamma))
            if lenC == 0:
                break
            else:
                if lenC > beta:
                    B = keepBetaBest(C,h, beta)
                else:
                    B = C
            lenlcs +=1
    stoptime = timeit.default_timer()
    return lenlcs, (stoptime - starttime)

def UBHH(name, _Beta):
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
                    lenlcs, tim = main(lenstrings[k], lenalphabet[j], fname, inalphabet, False, _Beta)
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
                    lenlcs, tim = main(lenstrings[k], lenalphabet[j], fname, inalphabet, False, _Beta)
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
                    lenlcs, tim  = main(lenstrings[k], lenalphabet[j], fname, inalphabet, False, _Beta)
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
                    x, tim = main(lenstrings[j], lenalphabet[i], fname, inalphabet, False, _Beta)
                    answers.append(x)
                    _timav.append(tim)
                print(str(lenalphabet[i]) + "_" + str(lenstrings[j]) + "_1000.het0.1", np.mean(answers), np.mean(_timav))

    if name == "SARS-CoV-2":
        inalphabet = ['A', 'C', 'G', 'T']
        lenalphabet = len(inalphabet)
        lenstrings = [10,20,30,40,50,60,70,80,90,100,110]
        for i in range(len(lenstrings)):
            fname = 'cov_' + str(lenstrings[i]) + "_" + str(lenalphabet)
            lenlcs, tim = main(lenstrings[k], inalphabet, fname, inalphabet, False, _Beta)
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
                    lenlcs, tim = main(lenstrings[j], lenalphabet[i], fname, inalphabet, False, _Beta)
                    answers.append(lenlcs)
                    _timav.append(tim)
                print("ES_" + str(lenstrings[j]) + "_" + str(lenalphabet[i]), np.mean(answers), np.mean(_timav))



if __name__ == "__main__":
    datasetName = sys.argv[1]
    beamWidth = int(sys.argv[2])
    UBHH(datasetName, beamWidth)