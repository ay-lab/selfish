params.f1 = ""
params.f2 = ""
params.o = ""
params.r = 0
params.b = ""
params.sz = 1.6
params.i = 10
params.d = ""
params.v = ""
params.lm = ""

b_arg = params.b == "" ? "" : "-b $params.b "
sz_arg = params.sz == 1.6 ? "" : "-sz $params.sz "
i_arg = params.i == 10 ? "" : "-i $params.i "
d_arg = params.d == "" ? "" : "-d $params.d "
v_arg = params.v == "" ? "" : "-v $params.v "
lm_arg = params.lm == "" ? "" : "-lm $params.lm "

contact_map_one = Channel.fromPath( params.f1 )
contact_map_two = Channel.fromPath( params.f2 )

process blastThemAll {
  input:
  file cm1 from contact_map_one
  file cm2 from contact_map_two
  output:
  file 'selfish.npy'
  stdout res

  """
    /selfish/selfish/selfish.py -f1 $cm1 -f2 $cm2 -o ./ -r $params.r $b_arg$sz_arg$i_arg$d_arg$v_arg$lm_arg
  """

}
res.println {it}