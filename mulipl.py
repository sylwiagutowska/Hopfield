r3=["A_{lm_3}^{\kk *}u_{l_3}^*(r)","B_{lm_3}^{\kk^*}\dot{u}_{l_3}^*(r)",\
    "C_{lm_3}^{\kk^*}R^{lo*}_{lm_3}(r)","D_{lm_3}^{\kk^*}\dot{R}^{lo*}_{lm_3}(r)"]
r4=["A_{lm_4}^{\kkp}u_{l_4}(r)","B_{lm_4}^{\kkp}\dot{u}_{l_4}(r)",\
    "C_{lm_4}^{\kkp}R^{lo}_{lm_4}(r)","D_{lm_4}^{\kkp}\dot{R}^{lo}_{lm_4}(r)"]
r5=["A_{lm_5}^{\kkp *}u_{l_5}^*(r')","B_{lm_5}^{\kkp *}\dot{u}_{l_5}^*(r')",\
    "C_{lm_5}^{\kkp^*}R^{lo*}_{lm_5}(r')","D_{lm_5}^{\kkp^*}\dot{R}^{lo*}_{lm_5}(r')"]
r6=["A_{lm_6}^{\kk }u_{l_6}(r')","B_{lm_6}^{\kk}\dot{u}_{l_6}(r')",\
    "C_{lm_6}^{\kk}R^{lo}_{lm_6}(r')","D_{lm_6}^{\kk}\dot{R}^{lo}_{lm_6}(r')"]
nr=0
r,rp="",""
for i in r3:
 for j in r4:
  for k in r5:
   for l in r6:
    w=i+j+k+l+'+'
    if r3[0] in w and r4[0] in w: r+= "r_int34[0],"
    if r3[0] in w and r4[1] in w: r+= "r_int34[1],"
    if r3[0] in w and r4[2] in w: r+= "r_int34[2],"
    if r3[0] in w and r4[3] in w: r+= "r_int34[3],"
    if r3[1] in w and r4[0] in w: r+= "r_int34[4],"
    if r3[1] in w and r4[1] in w: r+= "r_int34[5],"
    if r3[1] in w and r4[2] in w: r+= "r_int34[6],"
    if r3[1] in w and r4[3] in w: r+= "r_int34[7],"
    if r3[2] in w and r4[0] in w: r+= "r_int34[8],"
    if r3[2] in w and r4[1] in w: r+= "r_int34[9],"
    if r3[2] in w and r4[2] in w: r+= "r_int34[10],"
    if r3[2] in w and r4[3] in w: r+= "r_int34[11],"
    if r3[3] in w and r4[0] in w: r+= "r_int34[12],"
    if r3[3] in w and r4[1] in w: r+= "r_int34[13],"
    if r3[3] in w and r4[2] in w: r+= "r_int34[14],"
    if r3[3] in w and r4[3] in w: r+= "r_int34[15],"
    nr+=1
    if nr%3==0: 
      w+="\\\ &"
      r+="\\\n"
print r
for i in r3:
 for j in r4:
  for k in r5:
   for l in r6:
    w=i+j+k+l+'+'
    if r5[0] in w and r6[0] in w: rp+= "r_int56[0],"
    if r5[0] in w and r6[1] in w: rp+= "r_int56[1],"
    if r5[0] in w and r6[2] in w: rp+= "r_int56[2],"
    if r5[0] in w and r6[3] in w: rp+= "r_int56[3],"
    if r5[1] in w and r6[0] in w: rp+= "r_int56[4],"
    if r5[1] in w and r6[1] in w: rp+= "r_int56[5],"
    if r5[1] in w and r6[2] in w: rp+= "r_int56[6],"
    if r5[1] in w and r6[3] in w: rp+= "r_int56[7],"
    if r5[2] in w and r6[0] in w: rp+= "r_int56[8],"
    if r5[2] in w and r6[1] in w: rp+= "r_int56[9],"
    if r5[2] in w and r6[2] in w: rp+= "r_int56[10],"
    if r5[2] in w and r6[3] in w: rp+= "r_int56[11],"
    if r5[3] in w and r6[0] in w: rp+= "r_int56[12],"
    if r5[3] in w and r6[1] in w: rp+= "r_int56[13],"
    if r5[3] in w and r6[2] in w: rp+= "r_int56[14],"
    if r5[3] in w and r6[3] in w: rp+= "r_int56[15],"
    nr+=1
    if nr%3==0: 
      w+="\\\ &"
      rp+="\\\n"
print rp

kk=""
for i in r3:
 for j in r4:
  for k in r5:
   for l in r6:
    w=i+j+k+l+'+'
    if r3[0] in w and r6[0] in w: kk+= "k_int36[0],"
    if r3[0] in w and r6[1] in w: kk+= "k_int36[1],"
    if r3[0] in w and r6[2] in w: kk+= "k_int36[2],"
    if r3[0] in w and r6[3] in w: kk+= "k_int36[3],"
    if r3[1] in w and r6[0] in w: kk+= "k_int36[4],"
    if r3[1] in w and r6[1] in w: kk+= "k_int36[5],"
    if r3[1] in w and r6[2] in w: kk+= "k_int36[6],"
    if r3[1] in w and r6[3] in w: kk+= "k_int36[7],"
    if r3[2] in w and r6[0] in w: kk+= "k_int36[8],"
    if r3[2] in w and r6[1] in w: kk+= "k_int36[9],"
    if r3[2] in w and r6[2] in w: kk+= "k_int36[10],"
    if r3[2] in w and r6[3] in w: kk+= "k_int36[11],"
    if r3[3] in w and r6[0] in w: kk+= "k_int36[12],"
    if r3[3] in w and r6[1] in w: kk+= "k_int36[13],"
    if r3[3] in w and r6[2] in w: kk+= "k_int36[14],"
    if r3[3] in w and r6[3] in w: kk+= "k_int36[15],"
    nr+=1
    if nr%3==0: 
      w+="\\\ &"
      kk+="\\\n"
print kk

kp=""
for i in r3:
 for j in r4:
  for k in r5:
   for l in r6:
    w=i+j+k+l+'+'
    if r4[0] in w and r5[0] in w: kp+= "k_int45[0],"
    if r4[0] in w and r5[1] in w: kp+= "k_int45[1],"
    if r4[0] in w and r5[2] in w: kp+= "k_int45[2],"
    if r4[0] in w and r5[3] in w: kp+= "k_int45[3],"
    if r4[1] in w and r5[0] in w: kp+= "k_int45[4],"
    if r4[1] in w and r5[1] in w: kp+= "k_int45[5],"
    if r4[1] in w and r5[2] in w: kp+= "k_int45[6],"
    if r4[1] in w and r5[3] in w: kp+= "k_int45[7],"
    if r4[2] in w and r5[0] in w: kp+= "k_int45[8],"
    if r4[2] in w and r5[1] in w: kp+= "k_int45[9],"
    if r4[2] in w and r5[2] in w: kp+= "k_int45[10],"
    if r4[2] in w and r5[3] in w: kp+= "k_int45[11],"
    if r4[3] in w and r5[0] in w: kp+= "k_int45[12],"
    if r4[3] in w and r5[1] in w: kp+= "k_int45[13],"
    if r4[3] in w and r5[2] in w: kp+= "k_int45[14],"
    if r4[3] in w and r5[3] in w: kp+= "k_int45[15],"
    nr+=1
    if nr%3==0: 
      w+="\\\ &"
      kp+="\\\n"
print kp
