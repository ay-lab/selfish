params.f1 = ""
params.f2 = ""
params.r = 0
params.b1 = "none1"
params.b2 = "none2"
params.sz = 1.6
params.i = 10
params.d = ""
params.v = ""
params.lm = ""
params.ch = "n"

bias_file1 = file(params.b1)
bias_file2 = file(params.b2)

sz_arg = params.sz == 1.6 ? "" : "-sz $params.sz "
i_arg = params.i == 10 ? "" : "-i $params.i "
d_arg = params.d == "" ? "" : "-d $params.d "
v_arg = params.v == "" ? "" : "-v $params.v "
lm_arg = params.lm == "" ? "" : "-lm $params.lm "
ch_arg = params.ch == "n" ? "" : "-ch $params.ch "

contact_map_one = Channel.fromPath( params.f1 )
contact_map_two = Channel.fromPath( params.f2 )

process run_selfish {
  input:
  file cm1 from contact_map_one
  file cm2 from contact_map_two
  file b_f1 from  bias_file1
  file b_f2 from  bias_file2
  output:
  file 'selfish.npy'
  stdout res

  script:
  def b_arg1 = b_f1.name != "none1" ? "-b1 $b_f1 " : ""
  def b_arg2 = b_f2.name != "none2" ? "-b2 $b_f2 " : ""
  """
    /selfish/selfish/selfish.py -f1 $cm1 -f2 $cm2 -o ./selfish.npy -r $params.r $b_arg1$b_arg2$ch_arg$sz_arg$i_arg$d_arg$v_arg$lm_arg
  """

}
res.println {it}
