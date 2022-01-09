function xi=PtPoints(collo1,collo2,collo3)
xi=zeros(81,4);
index1=1;
index2=2;
stop=3;
int_counter=1;
counter=1;
repeat=0;
col=1;
while repeat <= 3
 if repeat~=0
 index1=index1*3;
 index2=index2*3;
 stop=stop*3;
 int_counter=1;
 counter=1;
 col=1+repeat;
 end
 while counter<=81
 if int_counter<=index1
 xi(counter,col)=collo1;
 int_counter=int_counter+1;
 counter=counter+1;
 elseif index1<=int_counter && int_counter<=index2
 xi(counter,col)=collo2;
 int_counter=int_counter+1;
 counter=counter+1;
 elseif index2<int_counter && int_counter<=stop
 xi(counter,col)=collo3;
 int_counter=int_counter+1;
 counter=counter+1;
 else
 int_counter=1;
 end
 end
 repeat=repeat+1;
end
end