function normkey(s, t){
  t=s
  gsub(/_M[0-9]+_B[0-9]+/,"",t)
  gsub(/_MB[0-9]+/,"",t)
  gsub(/_B[0-9]+/,"",t)
  gsub(/_(launch|run[0-9]+|clean|c[0-9]+|pilot|queue|onecore|debug|new)$/,"",t)
  gsub(/__+/,"_",t)
  sub(/_$/,"",t)
  return t
}
{
  name=$0; M=-1; B=-1;
  if (match(name,/_M([0-9]+)_B([0-9]+)/,a)) {M=a[1]+0; B=a[2]+0}
  else if (match(name,/_MB([0-9]+)/,a)) {M=a[1]+0; B=a[1]+0}
  else if (match(name,/_B([0-9]+)/,a)) {M=0; B=a[1]+0}
  else next
  key=normkey(name)
  qual=((M>=500 && B>=500)?1:0)
  print key"\t"name"\t"M"\t"B"\t"qual
}
