T = 2.5 ;   J = 1 ;  K = 1 ; B = 0 ;
N = 1e7  ;  N_equil = 1e6 ; Dim = [8 8];

is_equil = true ; is_biased = true ;

row = 1 ; col = 2 ;

e_hist = -J*(2*Dim(row)*Dim(col)) : 4*J : J*(2*Dim(row)*Dim(col)) ;

h = waitbar(0,'Please wait...');


m_equil = zeros(1,N) ;
e_equil = zeros(1,N) ;

m_hist = -Dim(row)*Dim(col) : 2 : Dim(row)*Dim(col) ;


if is_equil == true
    init_conf =  ones(Dim(row),Dim(col))  ;
    e =  -J*(2*Dim(row)*Dim(col)) ;
    
    is_biased = false ;
end

if is_biased == false
    
    e_tmn = [e_hist' zeros(length(e_hist),5) ] ;
    e_tm = [e_hist' zeros(length(e_hist),5) ] ;
    
    
end

m_hist = [ m_hist'  zeros(length(m_hist),1) ] ;
e_hist = [ e_hist'  zeros(length(e_hist),1) ]; 

conf = init_conf ;

n_acc = 0 ;



for iter_no = 1 : N
    
    spin = [ randi(Dim(row)), randi(Dim(col)) ] ;
    old_spin = conf(spin(row),spin(col)) ; e_old = e ;
    new_spin = - old_spin ;
    
    
    above = mod(spin(row) - 1 - 1, size(conf,1)) + 1; above = conf(above,spin(col)) ;
    below = mod(spin(row) + 1 - 1, size(conf,1)) + 1; below = conf(below,spin(col)) ;
    left  = mod(spin(col) - 1 - 1, size(conf,2)) + 1; left = conf(spin(row),left) ;
    right = mod(spin(col) + 1 - 1, size(conf,2)) + 1; right = conf(spin(row),right) ;
    
    neigh = [above below left right] ;
    
    % dE = J * ( new_spin - old_spin ) * ((above + below + left + right) - B ) ;
    bias = 0 ;             e_old_index = find(e_tm(:,1) == e_old,1) ;
    
    dE = J * (  - new_spin + old_spin ) * (sum(neigh) - B ) ;
    e_new = e + dE ;
    if is_biased == true && iter_no > 1000000
        if e_old - e_new == -4
            if isnan(e_tmn(e_old_index + 1,5)) == false
                bias = log(e_tmn(e_old_index + 1,5)) - log(e_tmn(e_old_index,3)) ;
            
            end
            
        elseif e_old - e_new == 4
            if isnan(e_tmn(e_old_index - 1,3)) == false
                bias = log(e_tmn(e_old_index - 1,3)) - log(e_tmn(e_old_index,5)) ;
            
            end
        elseif e_old - e_new == -8
            if isnan(e_tmn(e_old_index + 2,6)) == false
                bias = log(e_tmn(e_old_index + 2,6)) - log(e_tmn(e_old_index,4)) ;
            
            end
            
        elseif e_old - e_new == 8
            if isnan(e_tmn(e_old_index - 2,4)) == false
                bias = log(e_tmn(e_old_index - 2,4)) - log(e_tmn(e_old_index,6)) ;
            
            end
            
        end
    end
    
    prob_log = - dE/(K*T) ; prob = min( exp(prob_log) , 1 ) ;
    
    bias_prob_log = bias + prob_log ; 
    
    acc_prob = rand() ; acc_prob_log = log(acc_prob) ;
    e_new = e + dE ;
    
    if acc_prob_log <=  bias_prob_log
        conf(spin(row),spin(col)) = new_spin ;
        n_acc = n_acc + 1 ;
        e = e_new ;
        
    end
    
    
    e_index = find(e_hist(:,1) == e,1);
   
    
    m = sum(sum(conf)) ;
    m_index = find(m_hist(:,1) == m);
   
    
    m_hist(m_index, 2) = m_hist(m_index, 2) + 1 ;
    
    e_hist(e_index, 2) = e_hist(e_index, 2) + 1 ;
   
    
    e_tm(e_old_index, 2) = e_tm(e_old_index, 2) + 1 ;
    if e_new - e_old  == 4
        e_tm(e_old_index, 3) = e_tm(e_old_index, 3) + prob ;
    elseif e_new - e_old  == 8
        e_tm(e_old_index, 4) = e_tm(e_old_index, 4) + prob ;
    elseif e_new - e_old  == -4
        e_tm(e_old_index, 5) = e_tm(e_old_index, 5) + prob ;
    elseif e_new - e_old  == -8
        e_tm(e_old_index, 6) = e_tm(e_old_index, 6) + prob ;
    end
    
    
    if iter_no > 100000 && mod(iter_no,100000) == 0
        e_tmn = e_tm ;
        e_tmn(:,3:6) = e_tmn(:,3:6)./e_tm(:,2) ;
    end
    
    if mod(iter_no,1000000) == 0
        
        waitbar(iter_no / N,h)
    end
end
init_conf = conf ;
close(h)
% plot(m_hist(:,1)/(Dim(row)*Dim(col)),log(m_hist(:,2)./m_hist(1,2)))
% plot(e_hist(:,1)/(Dim(row)*Dim(col)),log(e_hist(:,2)./e_hist(1,2)))