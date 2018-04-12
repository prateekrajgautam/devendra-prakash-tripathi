clc
clear all
close all
%% user defined variables starts
format shorte
hold on
plotcolor=['-b^';'-g*';'-ro';'-cx';'-ks';'-yd';'-mp';'-wh';'-b^';'-g*';'-ro';'-cx';'-ks';'-yd';'-mp';'-wh'];
% datafrom=1;% 1 it will work fast 2 it will be more precise but take more time
int=1;
for n1=[8,16,32]
    for it1=[5]
        for q=2    
        %     if datafrom==1
                block=4;
                n=n1                %no.of users
                m=512;              %data length
                sl=64;              %spread length
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
        %         data=randint(n,m,[1,0]);                %generation of random bipolar data    
                data=randi([0,1],n,m);
                chip=spreador(spreadingunipolar,data);   
                transmit1=interleavor(chip,scrambrule);       %transmitting data interleavor   
                %% bpsk encoder to transmit & apply energy profile
                t2=2*((transmit1)>0)-1;
                [r,M]=size(t2);
                N=q*M;
                for i=1:r
                    X(i,:)=fft(t2(i,:));
                end
                for i=1:r
                    Y(i,:)=zeros(1,N);
                    Y(i,[1:q:N])=X(i,:);
                end
                for i=1:r
                    y1(i,:)=ifft(Y(i,:));
                end
                tx2=sum(y1);
                R=channelnoise(tx2,ebno);%awgn channel
                Y=fft(R);
                X=Y(:,[1:q:N]);
                r1=ifft(X);
                appllr=idmadecoderbpsk(sigma,hk,n,m,sl,itnum,chiplen,r1,spreading,scrambrule,useldpc);
            %%                     reciever
                e=0;
                appllrf=appllr>0;% decision whether 1 or 0 is recieved
                decodeddata=appllrf;
                [e,bertemp]=errortx(data,decodeddata);% check error
                error=error+e;       
                ber=(error/(n*m*bloc));
        %         [z,bloc,ebno,ber]
                end
                ebn(1,z)=ebno;
                ber1(1,z)=ber;
        end
        result(int).ebn=ebn;
        result(int).ber=ber1;
           if q==1
               lname{int}=['IDMA,n=' num2str(n) ',it=' num2str(itnum)]
           end
           if q>1
               lname{int}=['IFDMA with q=' num2str(q) ',n=' num2str(n) ',it=' num2str(itnum)];
           end
           int=int+1;
            save('temp1','q','n1','it1','int','result','plotcolor','lname')
                clear all
                load temp1
        end     
    end
end
    clc
    figure
    hold on
    for i=1:(int-1)
        semilogy(result(i).ebn,result(i).ber,plotcolor(i,:))     
        disp(lname{i})
        disp(result(i))
    end
     xlabel('Eb/No')
    ylabel('Bit Error Rate')
    grid on
    title('BER of distributed IFDMA for different no of users')
    legend(lname)
    saveas(gcf,'graph n','jpg')