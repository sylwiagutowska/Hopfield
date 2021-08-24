x=[["\\frac{\\partial V_{LM}}{\\partial r}","\\sin\\theta","\\cos\\phi Y_{LM}Y_3^*Y_4"],
   ["\\frac{V_{LM}}{r}","\\frac{\\partial Y_{LM}}{\\partial \\theta}\\cos\\theta","\\cos\\phi Y_3^*Y_4"],
   ["-\\frac{V_{LM}}{r}","\\frac{1}{\\sin\\theta}","\\frac{\\partial Y_{LM}}{\\partial \\phi}\\sin\\phi Y_3^*Y_4"]]


xp=[["\\frac{\\partial V'_{LM}}{\\partial r'}","\\sin\\theta'","\\cos\\phi' Y'_{LM}Y_5'^*Y_6"],
   ["\\frac{V'_{LM}}{r'}","\\frac{\\partial Y'_{LM}}{\\partial \\theta'}\\cos\\theta'","\\cos\\phi' Y_5'^*Y_6"],
   ["-\\frac{V'_{LM}}{r'}","\\frac{1}{\\sin\\theta'}","\\frac{\\partial Y'_{LM}}{\\partial \\phi'}\\sin\\phi' Y_5'^*Y_6"]]

y=[["\\frac{\\partial V_{LM}}{\\partial r}","\\sin\\theta","\\sin\\phi Y_{LM}Y_3^*Y_4"],
   ["\\frac{V_{LM}}{r}","\\frac{\\partial Y_{LM}}{\\partial \\theta}\\cos\\theta","\\sin\\phi Y_3^*Y_4"],
   ["\\frac{V_{LM}}{r}","\\frac{1}{\\sin\\theta}","\\frac{\\partial Y_{LM}}{\\partial \\phi}\\cos\\phi Y_3^*Y_4"]]

yp=[["\\frac{\\partial V'_{LM}}{\\partial r'}","\\sin\\theta'","\\sin\\phi'  Y'_{LM}Y_5'^*Y_6"],
   ["\\frac{V'_{LM}}{r'}","\\frac{\\partial  Y'_{LM}}{\\partial \\theta'}\\cos\\theta'","\\sin\\phi' Y_5'^*Y_6"],
   ["\\frac{V'_{LM}}{r'}","\\frac{1}{\\sin\\theta'}","\\frac{\\partial  Y'_{LM}}{\\partial \\phi'}\\cos\\phi' Y_5'^*Y_6"]]

z=[["\\frac{\\partial V_{LM}}{\\partial r}","\\cos\\theta Y_{LM}Y_3^*Y_4",""],
   ["-\\frac{V_{LM}}{r}","\\frac{\\partial Y_{LM}}{\\partial \\theta}\\sin\\theta Y_3^*Y_4",""],
     ["","",""] ]

zp=[ ["\\frac{\\partial V'_{LM}}{\\partial r'}","\\cos\\theta' Y'_{LM}Y_5'^*Y_6",""],
     ["-\\frac{V'_{LM}}{r'}","\\frac{\\partial Y'_{LM}}{\\partial \\theta'}\\sin\\theta'Y_5'^*Y_6",""],
     ["","",""] ]

dot=''

for i in range(3):
 for j in range(3):
  for k in range(3):
    dot+=x[i][k]+xp[j][k]
  dot+="+\n"
 dot+="\\\\ \n&\n"
print dot
print
dot=''
for i in range(3):
 for j in range(3):
  for k in range(3):
    dot+=y[i][k]+yp[j][k]
  dot+="+\n"
 dot+="\\\\ \n&\n"
print dot
print
dot=''
for i in range(2):
 for j in range(2):
  for k in range(3):
    dot+=z[i][k]+zp[j][k]
  dot+="+\n"
 dot+="\\\\ \n&\n"
#print dot.replace("+-","-")


