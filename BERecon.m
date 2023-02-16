
a=4000;
durNum = 10
norms = zeros(durNum, 1);
load forcings.txt
for l = 1:durNum
    
    load B1.txt;
    g = 1; 
    V_Re = zeros(2000, 1);
    V_T = ones(2000, 1);
    B = B1;
    [height,width] = size(B);
    
    % the value 0.2 is the m_0 normalization from within the init_vals in
  
    % no equal to P(Bmatrixconnection)*image_vec_len
    for i=1:(height/2)
        for j=1:width
            B(i,j) = B(i,j) * 1.25 * 1/152  * 0.2;
        end
    end
    
    for i=height/2:height
        for j=1:width
            B(i,j) = B(i,j) * 1.0 * 1/152  * 0.2;
        end
    end
    
   
    mtx = dctmtx(width^.5);
    k =kron(mtx', mtx');
    
  
    N = B*k;

    % load f_rate_1.txt;
    f_str = ['f_rate_' num2str(l) '.txt'];
    fid1 = fopen(f_str, 'r');
		f_rate = textscan(fid1, '%f', a);
		fclose(fid1);
        m = f_rate{1}(1:2000);
        m2 = f_rate{1}(2001:4000);

    % load avg_volages_1.txt;
    v_str = ['avg_voltages_' num2str(l) '.txt'];
    fid2 = fopen(v_str, 'r');
		avg_voltages = textscan(fid2, '%f', a);
		fclose(fid2);
        v = avg_voltages{1}(1:2000);

    % load avg_neur_thresholds_1.txt;
    t_str = ['avg_neur_thresholds_' num2str(l) '.txt'];
    fid3 = fopen(t_str, 'r');
		avg_neur = textscan(fid3, '%f', a);
		fclose(fid3);
        V_T = avg_neur{1}(1:2000);


    
    
    load connectivity.txt;
    R = connectivity(1:2000, 1:2000);
    R2 = connectivity(1:2000, 2001:4000);
    
    %solve L1 optimization problem
    b=-g*(V_Re-v) - R*m - R2*m2 + m.*(V_T - V_Re);    %RHS of input-output map
    
    x2 = (cs_omp(b,N,width))';
    
    
    %convert image vector back to pixel mtx
    S = size(x2);                     
    D = zeros(width^(.5),width^(.5));
    k=1;
    j=1;
    for i = 1 : width
        if(k > width^(.5))
            k=1;
            j=j+1;
        end
        D(k,j) = x2(i);
        k=k+1;
    end
    
    
    %Invert 2d discrete cosine transform
   
    D = mtx'*D*mtx;
    
    save('fig2.txt', 'D', '-ascii')
    
    

    
    imagesc(uint8(D));
    colormap(gray);
    axis off    %need these after!
    axis image 
    axis square
    set(gcf,'Position',[200,200,200,200]);
    
    load stripes.jpg
    im1 = imread('stripes.jpg');
    Im = im1;
    %convert to useful form

    if(ndims(Im)==3)
        im1 = rgb2gray(Im);
    else
        im1 = double(Im);
    end

    S = size(im1);                     
    IM = zeros(width^(.5),width^(.5));
    k=1;
    j=1;
    for i = 1 : width
        if(k > width^(.5))
            k=1;
            j=j+1;
        end
        IM(k,j) = im1(i);
        k=k+1;
    end
    

    
 
    norms(l, 1) = norm(D - double(IM), 'fro')/norm(double(IM),'fro') 
end
