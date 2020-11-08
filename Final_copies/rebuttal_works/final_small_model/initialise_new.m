function [lissom] = initialise_new(lissom,opts,X1)  %initialize afferent and lateral weights(excitation and inhibition)
for ii = 2 : numel(lissom.layers)% for each layer
    clear aff_sheet
    if ii==2
        aff_sheet=zeros(size(X1(:,:,1)));
%     elseif ii==3
%         aff_sheet=zeros(lissom.layers{ii-1}.dim);
    end
    [lissom]=prefix_mats(lissom,aff_sheet,ii,opts);
    aff_sheet=lissom.layers{ii}.aff_sheet;
    WAFF=zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(1),numel(aff_sheet));
    for ii1=1:size(lissom.layers{ii}.I,1)
        waff_temp=rand(size(lissom.layers{ii}.I,2),1);% afferent weight initialization
        WAFF(ii1,lissom.layers{ii}.I(ii1,:))=waff_temp/sum(sum(abs(waff_temp)));
    end
    %     keyboard;
    lissom.layers{ii}.waff=sparse(WAFF);
    
    count=0;
    Wexc=zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(1));
    Winh=zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(1));
    %     Iexc=zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(1),(2*lissom.layers{ii}.exe_rad(1))^2);
    %     Iinh=zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(1),(2*lissom.layers{ii}.inhb_rad(1))^2-(2*lissom.layers{ii}.exe_rad(1))^2);
    I_exc_master=zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2));
    I_inh_master=zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2));
    [xgrid, ygrid] = meshgrid(1:lissom.layers{ii}.dim(2), 1:lissom.layers{ii}.dim(1));    % 10 X 10 (1 to 10) each
    if ii==3
        sigma(1)=0.5*lissom.layers{ii-1}.rf_fb(1);
        sigma(2)=0.5*lissom.layers{ii-1}.rf_fb(2);
        G=zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2));
    end
    for k = 1 : lissom.layers{ii}.dim(1)
        %each row in layer     %1 - 10
        for j = 1 : lissom.layers{ii}.dim(2)
            %              if k==8
            %             keyboard;
            %         end
            Iextemp=zeros(1,lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2));
            Iintemp=zeros(1,lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2));
            count=count+1;
            x = xgrid - k;    % offset the current position% Make a logical circular region set to 1, the rest to zero
            y = ygrid - j;    % 10 X 10 (1 to 10)
            nbrexc = x.^2 + y.^2 <= lissom.layers{ii}.exe_rad.^2;     %exci circle rad = 4
            outer_inhb = x.^2 + y.^2 <= lissom.layers{ii}.inhb_rad.^2;%inhi circle rad = 10
            nbrinhb= ~nbrexc.*outer_inhb;
            
            inhb_idx=(nbrinhb)>0;  %to find 1
            exc_idx=(nbrexc)>0;    %to find 1
            %             Iexc(count,1:numel(find((nbrexc)>0)))=find((nbrexc)>0);
            %             Iinh(count,1:numel(find((nbrinhb)>0)))=find((nbrinhb)>0);
            Iextemp((nbrexc)>0)=1;
            I_exc_master(count,:)=Iextemp;
            Iintemp((nbrinhb)>0)=1;
            I_inh_master(count,:)=Iintemp;
            
            x=rand(size(Iextemp(nbrexc>0),2),1);
            Wexc(count,exc_idx)=x;
            Wexc(count,:)= Wexc(count,:)/sum(sum(abs(Wexc(count,:))));
            x1=rand(size(Iintemp((nbrinhb)>0),2),1);
            Winh(count, inhb_idx)=x1;
            Winh(count,:)= Winh(count,:)/sum(sum(abs(Winh(count,:))));
            if ii==3
                xx3 = (((xgrid-k) + (ygrid-j))/sigma(1)).^2;
                yy3 = (((xgrid-k) - (ygrid-j))/sigma(2)).^2;
                gtemp = exp(-(xx3+yy3)/2);
                gtemp(gtemp<1e-6)=0;
                G(count,:)=reshape(gtemp,1,lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2));
                if opts.laton==0
                    Wexc(count,:)=0;
                    Winh(count,:)=0;
                    Wexc(count,count)=rand;
                end
            end
        end
    end
    
    %     lissom.layers{ii}.PreSynLat = zeros(lissom.layers{ii}.dim); %10 x 10 lateral response
    lissom.layers{ii}.Zprev =zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2),1);
    lissom.layers{ii}.Zold =zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2),1);
    lissom.layers{ii}.I_exc_master=sparse(I_exc_master);
    lissom.layers{ii}.I_inh_master=sparse(I_inh_master);
    lissom.layers{ii}.wexc=sparse(Wexc);
    lissom.layers{ii}.winh=sparse(Winh);
%     %       lissom.layers{ii}.Iexc=Iexc;
%     %       lissom.layers{ii}.Iinh=Iinh;
%     if ii==2
%         lissom.layers{ii}.threshold = 0.08 * ones(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2),1);
%     end
%     if ii==3
%         lissom.layers{ii}.Gaussian_mat=G;
%         lissom.layers{ii}.threshold = 0.01 * ones(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2),1);
%         
%     end
%     lissom.layers{ii}.rectify_res = zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2),1);
%     lissom.layers{ii}.ave_res = zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2),1);
%     lissom.layers{ii}.prev_ave_res = 0.024*ones(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2),1);
end
end