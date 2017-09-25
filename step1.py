# coding=utf-8
import cgi, sys, math, pp, numpy, random, copy

UnrPlinkPath = sys.argv[1]
ExrPlinkPath = sys.argv[2]
OutPlinkPath = sys.argv[3]
maxsites = int(sys.argv[4])
ncpus = int(sys.argv[5])

# ================== cubex ====================
# Readin: an 9-len list
# Return: gametes frq
# =============================================
def cubex(inputarray):
    result = {}
    n1111 = float(inputarray[0])
    n1112 = float(inputarray[1])
    n1122 = float(inputarray[2])
    n1211 = float(inputarray[3])
    n1212 = float(inputarray[4])
    n1222 = float(inputarray[5])
    n2211 = float(inputarray[6])
    n2212 = float(inputarray[7])
    n2222 = float(inputarray[8])
    n = (n1111 + n1112 + n1122 + n1211 + n1212 + n1222 + n2211 + n2212 + n2222)
    result["n"] = n
    p = ((n1111+n1112+n1122)*2.0+(n1211+n1212+n1222))/(2.0 * n)
    q = ((n1111+n1211+n2211)*2.0+(n1112+n1212+n2212))/(2.0 * n)
    n11 = (2.0*n1111 + n1112 + n1211) 
    n12 = (2.0*n1122 + n1112 + n1222)
    n21 = (2.0*n2211 + n2212 + n1211)
    n22 = (2.0*n2222 + n2212 + n1222) 
    a0 = -n11*p*q
    a1 = -n11*(1.0 - 2.0*p - 2.0*q) - n1212*(1.0 - p - q) + 2.0*n*p*q
    a2 = 2.0*n*(1.0 - 2.0*p - 2.0*q) - 2.0*n11 - n1212
    a3 = 4.0 * n
    minhap = n11 / (2.0 * float(n))
    maxhap = (n11 + n1212) / (2.0 * float(n))
    result["minhap"] = minhap
    result["maxhap"] = maxhap
    if p < 1.0 and p > 0.0:
        result["hwchisnp1"] = ((((n1111 + n1112 + n1122) - ((p ** 2)*n)))**2)/((p ** 2)*n)+\
                    ((((n1211 + n1212 + n1222) - ((2 * p * (1.0-p))*n)))**2)/((2 * p * (1.0-p))*n)+\
                    ((((n2211 + n2212 + n2222) - (((1-p) ** 2)*n)))**2)/(((1-p) ** 2)*n)
    else:
        result["hwchisnp1"] = 0.0
    if q < 1.0 and q > 0.0:
        result["hwchisnp2"] = ((((n1111 + n1211 + n2211) - ((q ** 2)*n)))**2)/((q ** 2)*n)+\
                    ((((n1112 + n1212 + n2212) - ((2 * q * (1.0-q))*n)))**2)/((2 * q * (1.0-q))*n)+\
                    ((((n1122 + n1222 + n2222) - (((1-q) ** 2)*n)))**2)/(((1-q) ** 2)*n)
    else:
        result["hwchisnp2"] = 0.0
    
    a = a3
    b = a2
    c = a1
    dee = a0
    
    xN = -b/(3.0*a)
    d2 = (math.pow(b,2)-3.0*a*c)/(9*math.pow(a,2))
    yN = a * math.pow(xN,3) + b * math.pow(xN,2) + c * xN + dee
    yN2 = math.pow(yN,2)

    h2 = 4 * math.pow(a,2) * math.pow(d2,3)
    result["realnoproblem"] = 0
    if abs(yN2-h2) <= 0.00001:
        result["realnoproblem"] = 1

    if yN2 > h2:
        # option 1
        number1 = 0.0
        number2 = 0.0
        if (1.0/(2.0*a)*(-yN + math.pow((yN2 - h2),0.5))) < 0:
            number1 = -math.pow(-(1.0/(2.0*a)*(-yN + math.pow((yN2 - h2),0.5))),1.0/3.0)
        else: number1 = math.pow((1.0/(2.0*a)*(-yN + math.pow((yN2 - h2),0.5))),1.0/3.0)
        if (1.0/(2.0*a)*(-yN - math.pow((yN2 - h2),0.5))) < 0:
            number2 = -math.pow(-(1.0/(2.0*a)*(-yN - math.pow((yN2 - h2),0.5))),1.0/3.0)
        else: number2 = math.pow((1.0/(2.0*a)*(-yN - math.pow((yN2 - h2),0.5))),1.0/3.0)
        alpha = xN + number1 + number2
        result["inputdata"] = str(inputarray)
        result["alpha"] = alpha
        result["beta"] = "Not a real root"
        result["gamma"] = "Not a real root"
        result["p"] = p
        result["q"] = q
        result["pq"] = str(p * q)
        if result["alpha"] >= minhap - 0.00001 and result["alpha"] <= maxhap + 0.00001:
            result["alphaposs"] = 1
            f11 = result["alpha"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["alphaDprime"] = round(Dprime,3)
            result["alpharsquared"] = round(rsquared,4)
            result["alphaf11"] = f11
            result["alphaf12"] = f12
            result["alphaf21"] = f21
            result["alphaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["alphaposs"] = 0
            else:
                result["alpha-e1111"] = n * f11**2
                result["alpha-e1112"] = 2 * n * f11 * f12
                result["alpha-e1122"] = n * f12**2
                result["alpha-e1211"] = 2 * n * f11 * f21
                result["alpha-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22 
                result["alpha-e1222"] = 2 * n * f12 * f22
                result["alpha-e2211"] = n * f21**2
                result["alpha-e2212"] = 2 * n * f21 * f22
                result["alpha-e2222"] = n * f22**2
                if result["alpha-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["alpha-e1111"])**2/result["alpha-e1111"])
                else: chisq1111 = 0.0
                if result["alpha-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["alpha-e1112"])**2/result["alpha-e1112"])
                else: chisq1112 = 0.0
                if result["alpha-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["alpha-e1122"])**2/result["alpha-e1122"])
                else: chisq1122 = 0.0
                if result["alpha-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["alpha-e1211"])**2/result["alpha-e1211"])
                else: chisq1211 = 0.0
                if result["alpha-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["alpha-e1212"])**2/result["alpha-e1212"])
                else: chisq1212 = 0.0
                if result["alpha-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["alpha-e1222"])**2/result["alpha-e1222"])
                else: chisq1222 = 0.0
                if result["alpha-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["alpha-e2211"])**2/result["alpha-e2211"])
                else: chisq2211 = 0.0
                if result["alpha-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["alpha-e2212"])**2/result["alpha-e2212"])
                else: chisq2212 = 0.0
                if result["alpha-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["alpha-e2222"])**2/result["alpha-e2222"])
                else: chisq2222 = 0.0
                result["alpha-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else:
            result["alphaposs"] = 0
        result["betaposs"] = 0
        result["gammaposs"] = 0

    elif yN2 == h2:
        # option 2
        delta = math.pow((yN/2.0*a),(1.0/3.0))
        result["inputdata"] = str(inputarray)
        result["alpha"] = xN + delta
        result["beta"] = xN + delta
        result["gamma"] = xN - 2.0*delta
        result["p"] = str(p)
        result["q"] = str(q)
        result["pq"] = str(p * q)
        if result["alpha"] >= minhap - 0.00001 and result["alpha"] <= maxhap + 0.00001:
            result["alphaposs"] = 1
            f11 = result["alpha"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["alphaDprime"] = round(Dprime,3)
            result["alpharsquared"] = round(rsquared,4)
            result["alphaf11"] = f11
            result["alphaf12"] = f12
            result["alphaf21"] = f21
            result["alphaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["alphaposs"] = 0
            else:
                result["alpha-e1111"] = n * f11**2
                result["alpha-e1112"] = 2 * n * f11 * f12
                result["alpha-e1122"] = n * f12**2
                result["alpha-e1211"] = 2 * n * f11 * f21
                result["alpha-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22 
                result["alpha-e1222"] = 2 * n * f12 * f22
                result["alpha-e2211"] = n * f21**2
                result["alpha-e2212"] = 2 * n * f21 * f22
                result["alpha-e2222"] = n * f22**2
                if result["alpha-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["alpha-e1111"])**2/result["alpha-e1111"])
                else: chisq1111 = 0.0
                if result["alpha-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["alpha-e1112"])**2/result["alpha-e1112"])
                else: chisq1112 = 0.0
                if result["alpha-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["alpha-e1122"])**2/result["alpha-e1122"])
                else: chisq1122 = 0.0
                if result["alpha-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["alpha-e1211"])**2/result["alpha-e1211"])
                else: chisq1211 = 0.0
                if result["alpha-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["alpha-e1212"])**2/result["alpha-e1212"])
                else: chisq1212 = 0.0
                if result["alpha-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["alpha-e1222"])**2/result["alpha-e1222"])
                else: chisq1222 = 0.0
                if result["alpha-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["alpha-e2211"])**2/result["alpha-e2211"])
                else: chisq2211 = 0.0
                if result["alpha-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["alpha-e2212"])**2/result["alpha-e2212"])
                else: chisq2212 = 0.0
                if result["alpha-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["alpha-e2222"])**2/result["alpha-e2222"])
                else: chisq2222 = 0.0
                result["alpha-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else:
            result["alphaposs"] = 0
        if result["beta"] >= minhap - 0.00001 and result["beta"] <= maxhap + 0.00001:
            result["betaposs"] = 1
            f11 = result["beta"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["betaDprime"] = round(Dprime,3)
            result["betarsquared"] = round(rsquared,4)
            result["betaf11"] = f11
            result["betaf12"] = f12
            result["betaf21"] = f21
            result["betaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["betaposs"] = 0
            else:
                result["beta-e1111"] = n * f11**2
                result["beta-e1112"] = 2 * n * f11 * f12
                result["beta-e1122"] = n * f12**2
                result["beta-e1211"] = 2 * n * f11 * f21
                result["beta-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22 
                result["beta-e1222"] = 2 * n * f12 * f22
                result["beta-e2211"] = n * f21**2
                result["beta-e2212"] = 2 * n * f21 * f22
                result["beta-e2222"] = n * f22**2
                if result["beta-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["beta-e1111"])**2/result["beta-e1111"])
                else: chisq1111 = 0.0
                if result["beta-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["beta-e1112"])**2/result["beta-e1112"])
                else: chisq1112 = 0.0
                if result["beta-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["beta-e1122"])**2/result["beta-e1122"])
                else: chisq1122 = 0.0
                if result["beta-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["beta-e1211"])**2/result["beta-e1211"])
                else: chisq1211 = 0.0
                if result["beta-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["beta-e1212"])**2/result["beta-e1212"])
                else: chisq1212 = 0.0
                if result["beta-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["beta-e1222"])**2/result["beta-e1222"])
                else: chisq1222 = 0.0
                if result["beta-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["beta-e2211"])**2/result["beta-e2211"])
                else: chisq2211 = 0.0
                if result["beta-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["beta-e2212"])**2/result["beta-e2212"])
                else: chisq2212 = 0.0
                if result["beta-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["beta-e2222"])**2/result["beta-e2222"])
                else: chisq2222 = 0.0
                result["beta-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else:
            result["betaposs"] = 0
        if result["gamma"] >= minhap - 0.00001 and result["gamma"] <= maxhap + 0.00001:
            result["gammaposs"] = 1
            f11 = result["gamma"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["gammaDprime"] = round(Dprime,3)
            result["gammarsquared"] = round(rsquared,4)
            result["gammaf11"] = f11
            result["gammaf12"] = f12
            result["gammaf21"] = f21
            result["gammaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["gammaposs"] = 0
            else:
                result["gamma-e1111"] = n * f11**2
                result["gamma-e1112"] = 2 * n * f11 * f12
                result["gamma-e1122"] = n * f12**2
                result["gamma-e1211"] = 2 * n * f11 * f21
                result["gamma-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22 
                result["gamma-e1222"] = 2 * n * f12 * f22
                result["gamma-e2211"] = n * f21**2
                result["gamma-e2212"] = 2 * n * f21 * f22
                result["gamma-e2222"] = n * f22**2
                if result["gamma-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["gamma-e1111"])**2/result["gamma-e1111"])
                else: chisq1111 = 0.0
                if result["gamma-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["gamma-e1112"])**2/result["gamma-e1112"])
                else: chisq1112 = 0.0
                if result["gamma-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["gamma-e1122"])**2/result["gamma-e1122"])
                else: chisq1122 = 0.0
                if result["gamma-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["gamma-e1211"])**2/result["gamma-e1211"])
                else: chisq1211 = 0.0
                if result["gamma-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["gamma-e1212"])**2/result["gamma-e1212"])
                else: chisq1212 = 0.0
                if result["gamma-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["gamma-e1222"])**2/result["gamma-e1222"])
                else: chisq1222 = 0.0
                if result["gamma-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["gamma-e2211"])**2/result["gamma-e2211"])
                else: chisq2211 = 0.0
                if result["gamma-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["gamma-e2212"])**2/result["gamma-e2212"])
                else: chisq2212 = 0.0
                if result["gamma-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["gamma-e2222"])**2/result["gamma-e2222"])
                else: chisq2222 = 0.0
                result["gamma-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else:
            result["gammaposs"] = 0
        
    elif yN2 < h2:
        #option 3
        h = math.pow(h2, 0.5)
        theta = ((math.acos(-yN/h))/3.0)
        #delta = math.pow((yN/2.0*a),(1.0/3.0)) # is it correct to reuse this?
        delta = math.pow(d2,0.5)
        result["inputdata"] = str(inputarray)
        result["alpha"] = xN + 2.0 * delta * math.cos(theta)
        result["beta"] = xN + 2.0 * delta * math.cos(2.0 * math.pi/3.0 + theta)
        result["gamma"] = xN + 2.0 * delta * math.cos(4.0 * math.pi/3.0 + theta)
        result["p"] = p
        result["q"] = q
        result["pq"] = p * q
        if result["alpha"] >= minhap - 0.00001 and result["alpha"] <= maxhap + 0.00001:
            result["alphaposs"] = 1
            f11 = result["alpha"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["alphaDprime"] = round(Dprime,3)
            result["alpharsquared"] = round(rsquared,4)
            result["alphaf11"] = f11
            result["alphaf12"] = f12
            result["alphaf21"] = f21
            result["alphaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["alphaposs"] = 0
            else:
                result["alpha-e1111"] = n * f11**2
                result["alpha-e1112"] = 2 * n * f11 * f12
                result["alpha-e1122"] = n * f12**2
                result["alpha-e1211"] = 2 * n * f11 * f21
                result["alpha-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22 
                result["alpha-e1222"] = 2 * n * f12 * f22
                result["alpha-e2211"] = n * f21**2
                result["alpha-e2212"] = 2 * n * f21 * f22
                result["alpha-e2222"] = n * f22**2
                if result["alpha-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["alpha-e1111"])**2/result["alpha-e1111"])
                else: chisq1111 = 0.0
                if result["alpha-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["alpha-e1112"])**2/result["alpha-e1112"])
                else: chisq1112 = 0.0
                if result["alpha-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["alpha-e1122"])**2/result["alpha-e1122"])
                else: chisq1122 = 0.0
                if result["alpha-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["alpha-e1211"])**2/result["alpha-e1211"])
                else: chisq1211 = 0.0
                if result["alpha-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["alpha-e1212"])**2/result["alpha-e1212"])
                else: chisq1212 = 0.0
                if result["alpha-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["alpha-e1222"])**2/result["alpha-e1222"])
                else: chisq1222 = 0.0
                if result["alpha-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["alpha-e2211"])**2/result["alpha-e2211"])
                else: chisq2211 = 0.0
                if result["alpha-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["alpha-e2212"])**2/result["alpha-e2212"])
                else: chisq2212 = 0.0
                if result["alpha-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["alpha-e2222"])**2/result["alpha-e2222"])
                else: chisq2222 = 0.0
                result["alpha-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else: result["alphaposs"] = 0
        if result["beta"] >= minhap - 0.00001 and result["beta"] <= maxhap + 0.00001:
            result["betaposs"] = 1
            f11 = result["beta"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["betaDprime"] = round(Dprime,3)
            result["betarsquared"] = round(rsquared,4)
            result["betaf11"] = f11
            result["betaf12"] = f12
            result["betaf21"] = f21
            result["betaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["betaposs"] = 0
            else:
                result["beta-e1111"] = n * f11**2
                result["beta-e1112"] = 2 * n * f11 * f12
                result["beta-e1122"] = n * f12**2
                result["beta-e1211"] = 2 * n * f11 * f21
                result["beta-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22
                result["beta-e1222"] = 2 * n * f12 * f22
                result["beta-e2211"] = n * f21**2
                result["beta-e2212"] = 2 * n * f21 * f22
                result["beta-e2222"] = n * f22**2
                if result["beta-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["beta-e1111"])**2/result["beta-e1111"])
                else: chisq1111 = 0.0
                if result["beta-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["beta-e1112"])**2/result["beta-e1112"])
                else: chisq1112 = 0.0
                if result["beta-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["beta-e1122"])**2/result["beta-e1122"])
                else: chisq1122 = 0.0
                if result["beta-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["beta-e1211"])**2/result["beta-e1211"])
                else: chisq1211 = 0.0
                if result["beta-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["beta-e1212"])**2/result["beta-e1212"])
                else: chisq1212 = 0.0
                if result["beta-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["beta-e1222"])**2/result["beta-e1222"])
                else: chisq1222 = 0.0
                if result["beta-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["beta-e2211"])**2/result["beta-e2211"])
                else: chisq2211 = 0.0
                if result["beta-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["beta-e2212"])**2/result["beta-e2212"])
                else: chisq2212 = 0.0
                if result["beta-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["beta-e2222"])**2/result["beta-e2222"])
                else: chisq2222 = 0.0
                result["beta-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else: result["betaposs"] = 0
        if result["gamma"] >= minhap - 0.00001 and result["gamma"] <= maxhap + 0.00001:
            result["gammaposs"] = 1
            f11 = result["gamma"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["gammaDprime"] = round(Dprime,3)
            result["gammarsquared"] = round(rsquared,4)
            result["gammaf11"] = f11
            result["gammaf12"] = f12
            result["gammaf21"] = f21
            result["gammaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["gammaposs"] = 0
            else:
                result["gamma-e1111"] = n * f11**2 
                result["gamma-e1112"] = 2 * n * f11 * f12
                result["gamma-e1122"] = n * f12**2
                result["gamma-e1211"] = 2 * n * f11 * f21
                result["gamma-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22 
                result["gamma-e1222"] = 2 * n * f12 * f22
                result["gamma-e2211"] = n * f21**2
                result["gamma-e2212"] = 2 * n * f21 * f22
                result["gamma-e2222"] = n * f22**2
                if result["gamma-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["gamma-e1111"])**2/result["gamma-e1111"])
                else: chisq1111 = 0.0
                if result["gamma-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["gamma-e1112"])**2/result["gamma-e1112"])
                else: chisq1112 = 0.0
                if result["gamma-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["gamma-e1122"])**2/result["gamma-e1122"])
                else: chisq1122 = 0.0
                if result["gamma-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["gamma-e1211"])**2/result["gamma-e1211"])
                else: chisq1211 = 0.0
                if result["gamma-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["gamma-e1212"])**2/result["gamma-e1212"])
                else: chisq1212 = 0.0
                if result["gamma-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["gamma-e1222"])**2/result["gamma-e1222"])
                else: chisq1222 = 0.0
                if result["gamma-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["gamma-e2211"])**2/result["gamma-e2211"])
                else: chisq2211 = 0.0
                if result["gamma-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["gamma-e2212"])**2/result["gamma-e2212"])
                else: chisq2212 = 0.0
                if result["gamma-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["gamma-e2222"])**2/result["gamma-e2222"])
                else: chisq2222 = 0.0
                result["gamma-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else: result["gammaposs"] = 0
        
    else: result["Error"] = "No answer"
    # Return the result as a dictionary array
    #
    if "alpha-chisq" in result and "beta-chisq" in result and "gamma-chisq" in result:
        chisqLst = [result["beta-chisq"], result["alpha-chisq"], result["gamma-chisq"]]
        x = chisqLst.index(min(chisqLst))
        if x==0:
            return [result["betaf11"],result["betaf12"],result["betaf21"],result["betaf22"],result["betaDprime"], result["betarsquared"]]
        elif x==1:
            return [result["alphaf11"],result["alphaf12"],result["alphaf21"],result["alphaf22"],result["alphaDprime"], result["alpharsquared"]]
        else:
            return [result["gammaf11"],result["gammaf12"],result["gammaf21"],result["gammaf22"],result["gammaDprime"], result["gammarsquared"]]
    #
    if "alpha-chisq" in result and "beta-chisq" in result:
        if result["alpha-chisq"] > result["beta-chisq"]:
            return [result["betaf11"],result["betaf12"],result["betaf21"],result["betaf22"],result["betaDprime"], result["betarsquared"]]
        else:
            return [result["alphaf11"],result["alphaf12"],result["alphaf21"],result["alphaf22"],result["alphaDprime"], result["alpharsquared"]]
    #
    if "beta-chisq" in result and "gamma-chisq" in result:
        if result["gamma-chisq"] > result["beta-chisq"]:
            return [result["betaf11"],result["betaf12"],result["betaf21"],result["betaf22"],result["betaDprime"], result["betarsquared"]]
        else:
            return [result["gammaf11"],result["gammaf12"],result["gammaf21"],result["gammaf22"],result["gammaDprime"], result["gammarsquared"]]
    #
    if "alpha-chisq" in result and "gamma-chisq" in result:
        if result["alpha-chisq"] > result["gamma-chisq"]:
            return [result["gammaf11"],result["gammaf12"],result["gammaf21"],result["gammaf22"],result["gammaDprime"], result["gammarsquared"]]
        else:
            return [result["alphaf11"],result["alphaf12"],result["alphaf21"],result["alphaf22"],result["alphaDprime"], result["alpharsquared"]]
    if "alpha-chisq" in result:
        return [result["alphaf11"],result["alphaf12"],result["alphaf21"],result["alphaf22"],result["alphaDprime"], result["alpharsquared"]]
    if "beta-chisq" in result:
        return [result["betaf11"],result["betaf12"],result["betaf21"],result["betaf22"],result["betaDprime"], result["betarsquared"]]
    if "gamma-chisq" in result:
        return [result["gammaf11"],result["gammaf12"],result["gammaf21"],result["gammaf22"],result["gammaDprime"], result["gammarsquared"]]

# ================== count81 ==================
# Readin: pair-ped and rsloc1, rsloc2, ref1, ref2
# Return: 81 kinds of count frq lst
# =============================================
def count81(pair_ped, rsloc1, rsloc2, ref1, ref2):
    singleman = [str(i)+str(j) for i in range(3) for j in range(3)]
    count_Ref = [i+j for i in singleman for j in singleman]
    count_Lst = [0 for i in count_Ref]
    pairLine = 0
    with open(pair_ped+'.ped', 'r') as pedfile:
        for lines in pedfile:
            if pairLine < 1:
                # geno1
                geno1 = lines.split("\t")[rsloc1+6]
                # geno2
                geno2 = lines.split("\t")[rsloc2+6]
                # 
                pergeno = str(2-geno1.count(ref1))+str(2-geno2.count(ref2))
                pairLine += 1
            else:
                # geno1
                geno1 = lines.split("\t")[rsloc1+6]
                # geno2
                geno2 = lines.split("\t")[rsloc2+6]
                # 
                pergeno += str(2-geno1.count(ref1))+str(2-geno2.count(ref2))
                count_Lst[count_Ref.index(pergeno)] += 1.
                pairLine = 0
    # print "\t".join(map(str, map(lambda x: x/sum(count_Lst), count_Lst)))
    return count_Lst

# ================== plink2cubex ==============
# Readin: Unr-ped and Unr-frq
# <args>: rsname 1 and 2
# ============================================= 
def plink2cubex(UnrPlinkPath, rsloc1, rsloc2, ref1, ref2):
    n1111 = 0.
    n1112 = 0.
    n1122 = 0.
    n1211 = 0.
    n1212 = 0.
    n1222 = 0.
    n2211 = 0.
    n2212 = 0.
    n2222 = 0.
    with open(UnrPlinkPath + '.ped', 'r') as pedfile:
        for lines in pedfile:
            # geno1
            geno1 = lines.split("\t")[rsloc1+6]
            # geno2
            geno2 = lines.split("\t")[rsloc2+6]
            # count
            if 2-geno1.count(ref1)==0 and 2-geno2.count(ref2)==0:
                n1111 += 1
            if 2-geno1.count(ref1)==0 and 2-geno2.count(ref2)==1:
                n1112 += 1
            if 2-geno1.count(ref1)==0 and 2-geno2.count(ref2)==2:
                n1122 += 1
            if 2-geno1.count(ref1)==1 and 2-geno2.count(ref2)==0:
                n1211 += 1
            if 2-geno1.count(ref1)==1 and 2-geno2.count(ref2)==1:
                n1212 += 1
            if 2-geno1.count(ref1)==1 and 2-geno2.count(ref2)==2:
                n1222 += 1
            if 2-geno1.count(ref1)==2 and 2-geno2.count(ref2)==0:
                n2211 += 1
            if 2-geno1.count(ref1)==2 and 2-geno2.count(ref2)==1:
                n2212 += 1
            if 2-geno1.count(ref1)==2 and 2-geno2.count(ref2)==2:
                n2222 += 1
        return [n1111, n1112, n1122, n1211, n1212, n1222, n2211, n2212, n2222]

# ================== D_i ======================
# Readin: gams-freq and r
# Return: derivation
# =============================================
def D_i(gametes_frq, r):
    g1 = gametes_frq[0]; g2=gametes_frq[1]; g3=gametes_frq[2]; g4=gametes_frq[3]
    d_i= [(-1.0*g1*g4 + 1.0*g2*g3)/(1.0*g1**2 + 1.0*g1*g2 + 1.0*g1*g3 + 2*g1*g4*(-0.5*r + 0.5) + 1.0*g2*g3*r),(1.0*g1*g4 - 1.0*g2*g3)/(1.0*g1*g2 + 1.0*g1*g4*r + 1.0*g2**2 + 2*g2*g3*(-0.5*r + 0.5) + 1.0*g2*g4),0,(1.0*g1*g4 - 1.0*g2*g3)/(1.0*g1*g3 + 1.0*g1*g4*r + 2*g2*g3*(-0.5*r + 0.5) + 1.0*g3**2 + 1.0*g3*g4),(-1.0*g1*g4 + 1.0*g2*g3)/(2*g1*g4*(-0.5*r + 0.5) + 1.0*g2*g3*r + 1.0*g2*g4 + 1.0*g3*g4 + 1.0*g4**2),0,0,0,0,(-0.5*g1*g4 + 0.5*g2*g3)/(0.5*g1**2 + 0.5*g1*g2 + 0.5*g1*g3 + 2*g1*g4*(-0.25*r + 0.25) + 0.5*g2*g3*r),0,(0.5*g1*g4 - 0.5*g2*g3)/(0.5*g1*g2 + 0.5*g1*g4*r + 0.5*g2**2 + 2*g2*g3*(-0.25*r + 0.25) + 0.5*g2*g4),(0.5*g1*g4 - 0.5*g2*g3)/(0.5*g1*g3 + 0.5*g1*g4*r + 2*g2*g3*(-0.25*r + 0.25) + 0.5*g3**2 + 0.5*g3*g4),0,(-0.5*g1*g4 + 0.5*g2*g3)/(2*g1*g4*(-0.25*r + 0.25) + 0.5*g2*g3*r + 0.5*g2*g4 + 0.5*g3*g4 + 0.5*g4**2),0,0,0,0,(-1.0*g1*g4 + 1.0*g2*g3)/(1.0*g1**2 + 1.0*g1*g2 + 1.0*g1*g3 + 2*g1*g4*(-0.5*r + 0.5) + 1.0*g2*g3*r),(1.0*g1*g4 - 1.0*g2*g3)/(1.0*g1*g2 + 1.0*g1*g4*r + 1.0*g2**2 + 2*g2*g3*(-0.5*r + 0.5) + 1.0*g2*g4),0,(1.0*g1*g4 - 1.0*g2*g3)/(1.0*g1*g3 + 1.0*g1*g4*r + 2*g2*g3*(-0.5*r + 0.5) + 1.0*g3**2 + 1.0*g3*g4),(-1.0*g1*g4 + 1.0*g2*g3)/(2*g1*g4*(-0.5*r + 0.5) + 1.0*g2*g3*r + 1.0*g2*g4 + 1.0*g3*g4 + 1.0*g4**2),0,0,0,(-0.5*g1*g4 + 0.5*g2*g3)/(0.5*g1**2 + 0.5*g1*g2 + 0.5*g1*g3 + 2*g1*g4*(-0.25*r + 0.25) + 0.5*g2*g3*r),(0.5*g1*g4 - 0.5*g2*g3)/(0.5*g1*g2 + 0.5*g1*g4*r + 0.5*g2**2 + 2*g2*g3*(-0.25*r + 0.25) + 0.5*g2*g4),0,0,0,0,(0.5*g1*g4 - 0.5*g2*g3)/(0.5*g1*g3 + 0.5*g1*g4*r + 2*g2*g3*(-0.25*r + 0.25) + 0.5*g3**2 + 0.5*g3*g4),(-0.5*g1*g4 + 0.5*g2*g3)/(2*g1*g4*(-0.25*r + 0.25) + 0.5*g2*g3*r + 0.5*g2*g4 + 0.5*g3*g4 + 0.5*g4**2),0,(2*g1*g4*(-0.5*g1**2 - 0.5*g1*g2 - 0.5*g1*g3 + 2*g1*g4*(0.5*r - 0.5) - 0.5*g2*g3*r + 1.0*g2*g3*(-0.5*r + 0.5)) + 2*g2*g3*(0.5*g1**2 + 0.5*g1*g2 + 0.5*g1*g3 - 0.5*g1*g4*r + 1.0*g1*g4*(-0.5*r + 0.5) + 1.0*g2*g3*r))/(2*g1*g4*(g1**2*(-0.5*r + 0.5) + 2*g1*g2*(-0.25*r + 0.25) + 2*g1*g3*(-0.25*r + 0.25) + 2*g1*g4*(-0.5*r + 0.5)**2 + 1.0*g2*g3*r*(-0.5*r + 0.5)) + 2*g2*g3*(0.5*g1**2*r + 0.5*g1*g2*r + 0.5*g1*g3*r + 1.0*g1*g4*r*(-0.5*r + 0.5) + 0.5*g2*g3*r**2)),(2*g1*g4*(0.5*g1**2 + 0.5*g1*g3 - 1.0*g1*g4*r + 2.0*g1*g4*(-0.5*r + 0.5) - 0.5*g2**2 + 2*g2*g3*(1.0*r - 0.5) - 0.5*g2*g4) + 2*g2*g3*(-0.5*g1**2 - 0.5*g1*g3 + 2*g1*g4*(1.0*r - 0.5) + 0.5*g2**2 - 1.0*g2*g3*r + 2.0*g2*g3*(-0.5*r + 0.5) + 0.5*g2*g4))/(2*g1*g4*(0.5*g1**2*r + 0.5*g1*g2 + 0.5*g1*g3*r + 2.0*g1*g4*r*(-0.5*r + 0.5) + g2**2*(-0.5*r + 0.5) + 2*g2*g3*(0.25*r**2 + (-0.5*r + 0.5)**2) + 2*g2*g4*(-0.25*r + 0.25)) + 2*g2*g3*(g1**2*(-0.5*r + 0.5) + 0.5*g1*g2 + 2*g1*g3*(-0.25*r + 0.25) + 2*g1*g4*(0.25*r**2 + (-0.5*r + 0.5)**2) + 0.5*g2**2*r + 2.0*g2*g3*r*(-0.5*r + 0.5) + 0.5*g2*g4*r)),(2*g1*g4*(0.5*g1*g2 + 1.0*g1*g4*r + 0.5*g2**2 - 0.5*g2*g3*r + 1.0*g2*g3*(-0.5*r + 0.5) + 0.5*g2*g4) + 2*g2*g3*(-0.5*g1*g2 - 0.5*g1*g4*r + 1.0*g1*g4*(-0.5*r + 0.5) - 0.5*g2**2 + 2*g2*g3*(0.5*r - 0.5) - 0.5*g2*g4))/(2*g1*g4*(0.5*g1*g2*r + 0.5*g1*g4*r**2 + 0.5*g2**2*r + 1.0*g2*g3*r*(-0.5*r + 0.5) + 0.5*g2*g4*r) + 2*g2*g3*(2*g1*g2*(-0.25*r + 0.25) + 1.0*g1*g4*r*(-0.5*r + 0.5) + g2**2*(-0.5*r + 0.5) + 2*g2*g3*(-0.5*r + 0.5)**2 + 2*g2*g4*(-0.25*r + 0.25))),(2*g1*g4*(0.5*g1**2 + 0.5*g1*g2 - 1.0*g1*g4*r + 2.0*g1*g4*(-0.5*r + 0.5) + 2*g2*g3*(1.0*r - 0.5) - 0.5*g3**2 - 0.5*g3*g4) + 2*g2*g3*(-0.5*g1**2 - 0.5*g1*g2 + 2*g1*g4*(1.0*r - 0.5) - 1.0*g2*g3*r + 2.0*g2*g3*(-0.5*r + 0.5) + 0.5*g3**2 + 0.5*g3*g4))/(2*g1*g4*(0.5*g1**2*r + 0.5*g1*g2*r + 0.5*g1*g3 + 2.0*g1*g4*r*(-0.5*r + 0.5) + 2*g2*g3*(0.25*r**2 + (-0.5*r + 0.5)**2) + g3**2*(-0.5*r + 0.5) + 2*g3*g4*(-0.25*r + 0.25)) + 2*g2*g3*(g1**2*(-0.5*r + 0.5) + 2*g1*g2*(-0.25*r + 0.25) + 0.5*g1*g3 + 2*g1*g4*(0.25*r**2 + (-0.5*r + 0.5)**2) + 2.0*g2*g3*r*(-0.5*r + 0.5) + 0.5*g3**2*r + 0.5*g3*g4*r)),(2*g1*g4*(-0.5*g1**2 + 2*g1*g4*(2.0*r - 1.0) + 0.5*g2**2 - 2.0*g2*g3*r + 4.0*g2*g3*(-0.5*r + 0.5) + 0.5*g3**2 - 0.5*g4**2) + 2*g2*g3*(0.5*g1**2 - 2.0*g1*g4*r + 4.0*g1*g4*(-0.5*r + 0.5) - 0.5*g2**2 + 2*g2*g3*(2.0*r - 1.0) - 0.5*g3**2 + 0.5*g4**2))/(2*g1*g4*(g1**2*(-0.5*r + 0.5) + 0.5*g1*g2 + 0.5*g1*g3 + 2*g1*g4*(0.5*r**2 + 2*(-0.5*r + 0.5)**2) + 0.5*g2**2*r + 4.0*g2*g3*r*(-0.5*r + 0.5) + 0.5*g2*g4 + 0.5*g3**2*r + 0.5*g3*g4 + g4**2*(-0.5*r + 0.5)) + 2*g2*g3*(0.5*g1**2*r + 0.5*g1*g2 + 0.5*g1*g3 + 4.0*g1*g4*r*(-0.5*r + 0.5) + g2**2*(-0.5*r + 0.5) + 2*g2*g3*(0.5*r**2 + 2*(-0.5*r + 0.5)**2) + 0.5*g2*g4 + g3**2*(-0.5*r + 0.5) + 0.5*g3*g4 + 0.5*g4**2*r)),(2*g1*g4*(-0.5*g1*g2 - 1.0*g1*g4*r + 2.0*g1*g4*(-0.5*r + 0.5) - 0.5*g2**2 + 2*g2*g3*(1.0*r - 0.5) + 0.5*g3*g4 + 0.5*g4**2) + 2*g2*g3*(0.5*g1*g2 + 2*g1*g4*(1.0*r - 0.5) + 0.5*g2**2 - 1.0*g2*g3*r + 2.0*g2*g3*(-0.5*r + 0.5) - 0.5*g3*g4 - 0.5*g4**2))/(2*g1*g4*(2*g1*g2*(-0.25*r + 0.25) + 2.0*g1*g4*r*(-0.5*r + 0.5) + g2**2*(-0.5*r + 0.5) + 2*g2*g3*(0.25*r**2 + (-0.5*r + 0.5)**2) + 0.5*g2*g4 + 0.5*g3*g4*r + 0.5*g4**2*r) + 2*g2*g3*(0.5*g1*g2*r + 2*g1*g4*(0.25*r**2 + (-0.5*r + 0.5)**2) + 0.5*g2**2*r + 2.0*g2*g3*r*(-0.5*r + 0.5) + 0.5*g2*g4 + 2*g3*g4*(-0.25*r + 0.25) + g4**2*(-0.5*r + 0.5))),(2*g1*g4*(0.5*g1*g3 + 1.0*g1*g4*r - 0.5*g2*g3*r + 1.0*g2*g3*(-0.5*r + 0.5) + 0.5*g3**2 + 0.5*g3*g4) + 2*g2*g3*(-0.5*g1*g3 - 0.5*g1*g4*r + 1.0*g1*g4*(-0.5*r + 0.5) + 2*g2*g3*(0.5*r - 0.5) - 0.5*g3**2 - 0.5*g3*g4))/(2*g1*g4*(0.5*g1*g3*r + 0.5*g1*g4*r**2 + 1.0*g2*g3*r*(-0.5*r + 0.5) + 0.5*g3**2*r + 0.5*g3*g4*r) + 2*g2*g3*(2*g1*g3*(-0.25*r + 0.25) + 1.0*g1*g4*r*(-0.5*r + 0.5) + 2*g2*g3*(-0.5*r + 0.5)**2 + g3**2*(-0.5*r + 0.5) + 2*g3*g4*(-0.25*r + 0.25))),(2*g1*g4*(-0.5*g1*g3 - 1.0*g1*g4*r + 2.0*g1*g4*(-0.5*r + 0.5) + 2*g2*g3*(1.0*r - 0.5) + 0.5*g2*g4 - 0.5*g3**2 + 0.5*g4**2) + 2*g2*g3*(0.5*g1*g3 + 2*g1*g4*(1.0*r - 0.5) - 1.0*g2*g3*r + 2.0*g2*g3*(-0.5*r + 0.5) - 0.5*g2*g4 + 0.5*g3**2 - 0.5*g4**2))/(2*g1*g4*(2*g1*g3*(-0.25*r + 0.25) + 2.0*g1*g4*r*(-0.5*r + 0.5) + 2*g2*g3*(0.25*r**2 + (-0.5*r + 0.5)**2) + 0.5*g2*g4*r + g3**2*(-0.5*r + 0.5) + 0.5*g3*g4 + 0.5*g4**2*r) + 2*g2*g3*(0.5*g1*g3*r + 2*g1*g4*(0.25*r**2 + (-0.5*r + 0.5)**2) + 2.0*g2*g3*r*(-0.5*r + 0.5) + 2*g2*g4*(-0.25*r + 0.25) + 0.5*g3**2*r + 0.5*g3*g4 + g4**2*(-0.5*r + 0.5))),(2*g1*g4*(2*g1*g4*(0.5*r - 0.5) - 0.5*g2*g3*r + 1.0*g2*g3*(-0.5*r + 0.5) - 0.5*g2*g4 - 0.5*g3*g4 - 0.5*g4**2) + 2*g2*g3*(-0.5*g1*g4*r + 1.0*g1*g4*(-0.5*r + 0.5) + 1.0*g2*g3*r + 0.5*g2*g4 + 0.5*g3*g4 + 0.5*g4**2))/(2*g1*g4*(2*g1*g4*(-0.5*r + 0.5)**2 + 1.0*g2*g3*r*(-0.5*r + 0.5) + 2*g2*g4*(-0.25*r + 0.25) + 2*g3*g4*(-0.25*r + 0.25) + g4**2*(-0.5*r + 0.5)) + 2*g2*g3*(1.0*g1*g4*r*(-0.5*r + 0.5) + 0.5*g2*g3*r**2 + 0.5*g2*g4*r + 0.5*g3*g4*r + 0.5*g4**2*r)),0,(-0.5*g1*g4 + 0.5*g2*g3)/(0.5*g1**2 + 0.5*g1*g2 + 0.5*g1*g3 + 2*g1*g4*(-0.25*r + 0.25) + 0.5*g2*g3*r),(0.5*g1*g4 - 0.5*g2*g3)/(0.5*g1*g2 + 0.5*g1*g4*r + 0.5*g2**2 + 2*g2*g3*(-0.25*r + 0.25) + 0.5*g2*g4),0,0,0,0,(0.5*g1*g4 - 0.5*g2*g3)/(0.5*g1*g3 + 0.5*g1*g4*r + 2*g2*g3*(-0.25*r + 0.25) + 0.5*g3**2 + 0.5*g3*g4),(-0.5*g1*g4 + 0.5*g2*g3)/(2*g1*g4*(-0.25*r + 0.25) + 0.5*g2*g3*r + 0.5*g2*g4 + 0.5*g3*g4 + 0.5*g4**2),0,0,0,(-1.0*g1*g4 + 1.0*g2*g3)/(1.0*g1**2 + 1.0*g1*g2 + 1.0*g1*g3 + 2*g1*g4*(-0.5*r + 0.5) + 1.0*g2*g3*r),(1.0*g1*g4 - 1.0*g2*g3)/(1.0*g1*g2 + 1.0*g1*g4*r + 1.0*g2**2 + 2*g2*g3*(-0.5*r + 0.5) + 1.0*g2*g4),0,(1.0*g1*g4 - 1.0*g2*g3)/(1.0*g1*g3 + 1.0*g1*g4*r + 2*g2*g3*(-0.5*r + 0.5) + 1.0*g3**2 + 1.0*g3*g4),(-1.0*g1*g4 + 1.0*g2*g3)/(2*g1*g4*(-0.5*r + 0.5) + 1.0*g2*g3*r + 1.0*g2*g4 + 1.0*g3*g4 + 1.0*g4**2),0,0,0,0,(-0.5*g1*g4 + 0.5*g2*g3)/(0.5*g1**2 + 0.5*g1*g2 + 0.5*g1*g3 + 2*g1*g4*(-0.25*r + 0.25) + 0.5*g2*g3*r),0,(0.5*g1*g4 - 0.5*g2*g3)/(0.5*g1*g2 + 0.5*g1*g4*r + 0.5*g2**2 + 2*g2*g3*(-0.25*r + 0.25) + 0.5*g2*g4),(0.5*g1*g4 - 0.5*g2*g3)/(0.5*g1*g3 + 0.5*g1*g4*r + 2*g2*g3*(-0.25*r + 0.25) + 0.5*g3**2 + 0.5*g3*g4),0,(-0.5*g1*g4 + 0.5*g2*g3)/(2*g1*g4*(-0.25*r + 0.25) + 0.5*g2*g3*r + 0.5*g2*g4 + 0.5*g3*g4 + 0.5*g4**2),0,0,0,0,(-1.0*g1*g4 + 1.0*g2*g3)/(1.0*g1**2 + 1.0*g1*g2 + 1.0*g1*g3 + 2*g1*g4*(-0.5*r + 0.5) + 1.0*g2*g3*r),(1.0*g1*g4 - 1.0*g2*g3)/(1.0*g1*g2 + 1.0*g1*g4*r + 1.0*g2**2 + 2*g2*g3*(-0.5*r + 0.5) + 1.0*g2*g4),0,(1.0*g1*g4 - 1.0*g2*g3)/(1.0*g1*g3 + 1.0*g1*g4*r + 2*g2*g3*(-0.5*r + 0.5) + 1.0*g3**2 + 1.0*g3*g4),(-1.0*g1*g4 + 1.0*g2*g3)/(2*g1*g4*(-0.5*r + 0.5) + 1.0*g2*g3*r + 1.0*g2*g4 + 1.0*g3*g4 + 1.0*g4**2)]
    return map(lambda x: round(x,5), d_i)

def solution(gams_frq, count_set):
    Lst = [sum(numpy.array(D_i(gams_frq, r/2000.))*count_set) for r in range(-1,1001)]
    rLst = [r/2000. for r in range(-1,1001)]
    if Lst[-1]*Lst[0] <= 0:
        for i in range(len(Lst)-1):
            if Lst[i]*Lst[i+1] < 0:
               res = -1.*Lst[i]*(rLst[i+1]-rLst[i])/(Lst[i+1]-Lst[i]) + rLst[i]
            elif Lst[i] == 0:
                res = rLst[i]
            elif Lst[i+1] == 0:
                res = rLst[i+1]
    elif Lst[0] <= 0:
        res = 'NA'
    elif Lst[-1] > 0:
        res = 0.5
    else:
        res = 0
    if res == 'NA':
        return 'NA'
    else:
        return round(res,8)

def CalLikelihood(UnrPlinkPath, ExrPlinkPath, maxsites, rsi, max_j):
    strs = ''
    rsj = rsi + maxsites
    if rsj <= max_j:
        refLst = []
        rsLst  = []
        rsLoc  = []
        with open(UnrPlinkPath + '.map', 'r') as frqfile:
            for line in frqfile:
                rsLoc.append(line[:-1].split('\t')[3])
        loc1= rsLoc[rsi]
        loc2= rsLoc[rsj]
        #
        with open(UnrPlinkPath + '.frq', 'r') as frqfile:
            for line in frqfile:
                refLst.append(line.split('\t')[1])
                rsLst.append(line.split('\t')[0])
        ref1= refLst[rsi]
        ref2= refLst[rsj]
        #
        inputarray= plink2cubex(UnrPlinkPath, rsi, rsj, ref1, ref2)
        cubex_value = cubex(inputarray)
        # gams_frq  = map(lambda x: round(x,5), cubex_value[0:4])
        gams_frq = cubex_value[0:4]
        if 0 in gams_frq:
            for i in range(len(gams_frq)):
                if gams_frq[i]==0:
                    gams_frq[i] = 0.000001
        count_set = numpy.array(count81(ExrPlinkPath, rsi, rsj, ref1, ref2))  # 返回整数
        # recombination estimation
        res_rate = solution(gams_frq, count_set)
        # turban
        Lst81 = map(lambda x: x/sum(count_set), count_set)
        sim_resSet = [res_rate]
        multinomLst= numpy.random.multinomial(sum(count_set), Lst81, 500)
        # print count_set, multinomLst[0]
        for i in range(len(multinomLst)):
            sim_countset = multinomLst[i]
            sim_resSet.append(solution(gams_frq, sim_countset))
        # for t in range(200):
        #     sim_countset = [0 for i in range(len(Lst81))]
        #     for i in range(int(sum(count_set))):
        #         pr = random.random()
        #         for j in range(len(Lst81)):
        #             if pr<=sum(Lst81[0:j+1]):
        #                 sim_countset[j] += 1
        #                 break 
        # ========================
        # ana and return a ture value
        # ========================
        D_val = abs(gams_frq[3]*gams_frq[0]-gams_frq[1]*gams_frq[2])
        R_val = cubex_value[-1]
        
        # Est_val = numpy.array([i for i in sim_resSet if i<=0.4]) # 去掉绝对问题点
        # Est_val = numpy.array([i for i in Est_val if i <= numpy.percentile(Est_val, 95) and  i >= numpy.percentile(Est_val, 5)]) # 
        # print len(dat), numpy.percentile(dat, 90) - numpy.percentile(dat, 10), D_val, R_val
        # if D_val >= 0.15 or R_val >= 0.6:
        #     if len(dat) >= 30 and numpy.percentile(dat, 90) - numpy.percentile(dat, 10) <= 0.04:
        #         est = numpy.median([i for i in dat if i >= numpy.percentile(dat, 10) and i<=numpy.percentile(dat, 90)])
        #         interval_len = numpy.percentile(dat, 90) - numpy.percentile(dat, 10)
        #     else:
        #         est = 'NaN'
        #         interval_len = 'NaN' # error 1 means valid num is few

        # elif D_val >= 0.125 or R_val >= 0.4:
        #     if len(dat) >= 80 and numpy.percentile(dat, 90) - numpy.percentile(dat, 10) <= 0.035:
        #         est = numpy.median([i for i in dat if i >= numpy.percentile(dat, 10) and i<=numpy.percentile(dat, 90)])
        #         interval_len = numpy.percentile(dat, 90) - numpy.percentile(dat, 10)
        #     else:
        #         est = 'NaN'
        #         interval_len = 'NaN' # error 2 means valid interval is few

        # else:
        #     if len(dat) >= 80 and numpy.percentile(dat, 90) - numpy.percentile(dat, 10) <= 0.01:
        #         est = numpy.median([i for i in dat if i >= numpy.percentile(dat, 10) and i<=numpy.percentile(dat, 90)])
        #         interval_len = numpy.percentile(dat, 90) - numpy.percentile(dat, 10)
        #     else:
        #         est = 'NaN'
        #         interval_len = 'NaN'
        # if D_val >= 0.15 or R_val >= 0.8:
        #     if len(Est_val)*1./len(sim_resSet) >= 0.2: # avoid cannot percentile
        #         upper = numpy.percentile(Est_val, 75)
        #         lower = numpy.percentile(Est_val, 25) 
        #         Est = numpy.median([i for i in Est_val if i<=upper and i>=lower])
        #         est = 25*math.log((1+2*Est)/(1-2*Est))
        #         interval_len = upper-lower
        #     else:
        #         est = 'NA'
        #         interval_len = 'NA'

        # else:
        #     if len(Est_val)*1./len(sim_resSet) >= 0.2: # avoid cannot percentile
        #         upper = numpy.percentile(Est_val, 75)
        #         lower = numpy.percentile(Est_val, 25)
        #         if upper - lower <= 0.005:
        #             Est = numpy.median([i for i in Est_val if i<=upper and i>=lower])
        #             est = 25*math.log((1+2*Est)/(1-2*Est))
        #             interval_len = upper-lower
        #         else:
        #             est = 'NA'
        #             interval_len = 'NA'
        #     else:
        #         est = 'NA'
        #         interval_len = 'NA'

        res =  [rsLst[rsi], rsLst[rsj], rsi, rsj, loc1, loc2, D_val, R_val, gams_frq[0], gams_frq[1], gams_frq[2], gams_frq[3]]
        res.extend(sim_resSet)
        strs += ','.join(map(str, res))
    return strs

refLst = []

with open(UnrPlinkPath + '.frq', 'r') as frqfile:
    for line in frqfile:
        refLst.append(line.split('\t')[1])

max_j = len(refLst)

OutputFile = open(OutPlinkPath, "w")

OutputFile.write('rs1,rs2,rank1,rank2,loc1,loc2,D,R2,g1,g2,g3,g4,'+','.join(['est'+str(i) for i in range(501)])+'\n')

ppservers = ()
if len(sys.argv) > 2:
    ncpus = int(sys.argv[2])
    # Creates jobserver with ncpus workers
    job_server = pp.Server(ncpus, ppservers=ppservers)
else:
    # Creates jobserver with automatically detected number of workers
    job_server = pp.Server(ppservers=ppservers)

# jobs=[job_server.submit(CalLikelihood, (UnrPlinkPath, ExrPlinkPath, maxsites, rsi, max_j), (plink2cubex, count81, cubex,  D_i, solution), ('cgi', 'math', 'numpy', 'random', 'copy')) for rsi in range(max_j-5, max_j-1)]
jobs=[job_server.submit(CalLikelihood, (UnrPlinkPath, ExrPlinkPath, maxsites, rsi, max_j), (plink2cubex, count81, cubex,  D_i, solution), ('cgi', 'math', 'numpy', 'random', 'copy')) for rsi in range(max_j)]

for job in jobs:
    x = job()
    # print x
    if x != '':
        OutputFile.write(x+'\n')

OutputFile.close()
