T = 2.5 ;   J = 1 ;  K = 1 ; B = 0 ;
N = 1e7  ;  N_equil = 1e6 ; Dim = [8 8];

is_equil = false; is_biased = true ;

row = 1 ; col = 2 ;

m_hist = -Dim(row)*Dim(col) : 2 : Dim(row)*Dim(col) ;

h = waitbar(0,'Please wait...');






if is_equil == true
    init_conf = sign( 0.5 - rand(Dim(row),Dim(col)) ) ;
    
    
    is_biased = false ;
end

if is_biased == false
    
    m_tmn = [m_hist' zeros(length(m_hist),1) zeros(length(m_hist),1) zeros(length(m_hist),1)] ;
    m_tm = [m_hist'  zeros(length(m_hist),1) zeros(length(m_hist),1) zeros(length(m_hist),1)] ;
    
    
end


m_hist = [ m_hist'  zeros(length(m_hist),1) ] ;

conf = init_conf ;

n_acc = 0 ;

m = sum(sum(conf)) ;

for iter_no = 1 : N
    
    spin = [ randi(Dim(row)), randi(Dim(col)) ] ;
    old_spin = conf(spin(row),spin(col)) ; m_old = m ;
    new_spin = - old_spin ;
    
    
    above = mod(spin(row) - 1 - 1, size(conf,1)) + 1; above = conf(above,spin(col)) ;
    below = mod(spin(row) + 1 - 1, size(conf,1)) + 1; below = conf(below,spin(col)) ;
    left  = mod(spin(col) - 1 - 1, size(conf,2)) + 1; left = conf(spin(row),left) ;
    right = mod(spin(col) + 1 - 1, size(conf,2)) + 1; right = conf(spin(row),right) ;
    
    neigh = [above below left right] ;
    
    % dE = J * ( new_spin - old_spin ) * ((above + below + left + right) - B ) ;
    bias = 0 ; m_old_index = find(m_tm(:,1) == m_old) ;
    if is_biased == true && iter_no > 1000000
        if old_spin < new_spin
            if isnan(m_tmn(m_old_index + 1,4)) == false
                bias = log(m_tmn(m_old_index + 1,4)) - log(m_tmn(m_old_index,3)) ;
            else
                bias = inf ;
            end
            
        else
            if isnan(m_tmn(m_old_index - 1,3)) == false
                bias = log(m_tmn(m_old_index - 1,3)) - log(m_tmn(m_old_index,4)) ;
            else
                bias = inf ;
            end
            
        end
    end
    
    dE = J * (  - new_spin + old_spin ) * (sum(neigh) - B ) ;
    prob_log = - dE/(K*T) ; prob = min( exp(prob_log) , 1 ) ;
    
    bias_prob_log = bias + prob_log ; 
    
    acc_prob = rand() ; acc_prob_log = log(acc_prob) ;
    
    
    if acc_prob_log <=  bias_prob_log
        conf(spin(row),spin(col)) = new_spin ;
        n_acc = n_acc + 1 ;
        
    end
    
    m = sum(sum(conf)) ;
    m_index = find(m_hist(:,1) == m);
   
    
    m_hist(m_index, 2) = m_hist(m_index, 2) + 1 ;
   
    
    m_tm(m_old_index, 2) = m_tm(m_old_index, 2) + 1 ;
    if old_spin < new_spin
        m_tm(m_old_index, 3) = m_tm(m_old_index, 3) + prob ;
    else
        m_tm(m_old_index, 4) = m_tm(m_old_index, 4) + prob ;
    end
    
    
    if iter_no > 100000 && mod(iter_no,100000) == 0
        m_tmn = m_tm ;
        m_tmn(:,3:4) = m_tmn(:,3:4)./m_tm(:,2) ;
    end
    
    if mod(iter_no,1000000) == 0
        
        waitbar(iter_no / N,h)
    end

end
init_conf = conf ;
close(h) ;

%     plot(m_hist(:,1),log(m_hist(:,2)./m_hist(1,2)))