bbdose = bbslowfunc(doses,injectiontimes)
bbdoseAlOH = ALOHfunc(doses,injectiontimes)
bolusbb = bbfunc(doses,injectiontimes)
safebb5 = bbfunc(safedose5,dailytimes)
safebb50 = bbfunc(safedose50,dailytimes)
bmbb= bbfunc(breastmilkdose,dailytimes)
fmbb= bbfunc(formulamilkdose,dailytimes)


# adjust everything for background
netbodyburden = backgroundbb+bbdose
netbodyburdenAlOH = backgroundbb+bbdoseAlOH
bolusbb = bolusbb + backgroundbb
safebb5 = safebb5 + backgroundbb
safebb50 = safebb50 + backgroundbb
bmbb= bmbb + backgroundbb
fmbb= fmbb + backgroundbb