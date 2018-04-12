clc
clear all
close all
%% user defined variables starts
format shorte
plotcolor=['-b^';'-g*';'-ro';'-cx';'-ks';'-yd';'-mp';'-wh';'-b^';'-g*';'-ro';'-cx';'-ks';'-yd';'-mp';'-wh'];
% datafrom=1;% 1 it will work fast 2 it will be more precise but take more time
int=1;
for type=[1,2]%1 for SC-FDMA-IDMA & 2 for OFDM-IDMA
    if type<2
        typename='SC-FDMA-IDMA'
    end
    if type>1
        typename='OFDM-IDMA'
    end
    for n1=[8]
        for it1=[7]
            for q=2    
                block=10;
                n=n1                %no.of users
                m=512;              %data length
                sl=16;              %spread length
                chiplen=m*sl;       %chip length
                itnum=it1;            %no of iteration
                ebnostart=0;        %step iteration
                ebnostep=1;
                ebnonum=10;
                useldpc=[];
                useenergyprofile=0;%set useenergyprofile to one(1) to use useenergyprofile  else set it to 2
                usemodem=1; %select modulation type '1 for psk 2 for dpsk'
                %% user defined variables ends
                hk=ones(1,n);%assume we know h_k required for idma decoder
                %% idma 
                spreadingunipolar=spreadingsequece(sl,[1,0]);%spreading sequence producing {1,0,1,0--------}
                spreading=2*spreadingunipolar-1;
                scrambrule=idmascramble(n,chiplen);  
                %% the simulation process begins
                for z=1:ebnonum
                    ebno=ebnostart+z*ebnostep;
                    snr(z)=(10.^(ebno/10))/sl;
                    sigma=sqrt(0.5/snr(z));
                    error=0;    
                    for bloc=1:block
                    %%                     transmitter
                        %% transmitter section begins
                        data=randi([0,1],n,m);
                        chip=spreador(spreadingunipolar,data);   
                        transmit1=interleavor(chip,scrambrule);       %transmitting data interleavor   
                        %% bpsk encoder to transmit & apply energy profile
                        t2=2*((transmit1)>0)-1;
                        %% idma multipleser (sum)
                        [r,M]=size(t2);
                        N=q*M;
                        for i=1:r
                            X(i,:)=fft(t2(i,:));
                        end
                        if type<2
                            for i=1:r
                                Y(i,:)=zeros(1,N);
                                Y(i,[1:q:N])=X(i,:);
                            end
                            for i=1:r
                                y1(i,:)=ifft(Y(i,:));
                            end
                            tx2=sum(y1);
                        end
                        if type>1
                            tx2=sum(X);
                        end
                        %% noise
                        R=channelnoise(tx2,ebno);%awgn channel
%                         R=tx2;%ideal channel
                        if type<2
                            Y=fft(R);
                            X=Y(:,[1:q:N]);
                            r1=ifft(X);
                        end
                        if type>1
                            r1=ifft(R);
                        end
                        appllr=idmadecoderbpsk(sigma,hk,n,m,sl,itnum,chiplen,r1,spreading,scrambrule,useldpc);
%%                      reciever
%%                      idma decoder
                        e=0;
                        appllrf=appllr>0;% decision whether 1 or 0 is recieved
                        decodeddata=appllrf;
                        [e,bertemp]=errortx(data,decodeddata);% check error
                        error=error+e;       
                        ber=(error/(n*m*bloc));
                        [' z=' num2str(z) ' blockno=' num2str(bloc) ' ebno=' num2str(ebno) ' ber=' num2str(ber)]
                    end
                    ebn(1,z)=ebno;
                    ber1(1,z)=ber;
                end
                result(int).ebn=ebn;
                result(int).ber=ber1;
                if q==1
                    lname{int}=[typename 'n=' num2str(n) ',it=' num2str(itnum) 'block=' num2str(bloc)]
                end           
                if q>1
                    lname{int}=[typename 'q=' num2str(q) ',n=' num2str(n) ',it=' num2str(itnum) 'block=' num2str(bloc)];
                end
                int=int+1;
                save('temp1','q','n1','it1','int','result','plotcolor','lname','type','bloc')
                clear all
                load temp1
            end     
        end
    end
end

save('saveresult')
%%  plotting result 
clc
close all
hold on
for i=1:(int-1)
    hold
    semilogy(result(i).ebn,(result(i).ber),plotcolor(i,:))     
    disp(lname{i})
    disp(result(i))
end
xlabel('Eb/No')
ylabel('Bit Error Rate')
grid on
title('COMPARISION OF SC-FDMA-IDMA & OFDM-IDMA')
legend(lname)
saveas(gcf,'graph COMPARISION OF SC-FDMA-IDMA & OFDM-IDMA','jpg')