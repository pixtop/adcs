function [V w res_ev, it_ev] = fortran_subspace_iter_ev(A, m, ver, percentage, eps, maxit)
    %
    % routine that computes a certain amount of eigenvalues
    % and eigenvectors of a matrix A using the subspace iteration method
    %
    % INPUT arguments:
    %
    % A          = the input matrix
    %
    % m          = subspace dimension
    %             
    % ver        = method
    %
    % percentage = the percentage of the trace to be retained
    %
    % eps        = threshold value to stop the subspace iteration
    %
    % maxit      = the maximum number of iterations
    %
    %
    %
    % OUTPUT arguments:
    %
    % V      = a matrix containing eigenvector
    %
    % w      = a vector containing the computed eigenvalues
    %
    % res_ev = a vector containing the norm of the residual for each eigenvalue 
    %
    % it_ev  = a vector containing the number of iterations for each eigenvalue
    %

    if((nargin ~= 6) || (nargout~=4))
        fprintf('Wrong number of input/output arguments\n\n\nUsage:          \n');  
        fprintf('INPUT arguments:                                            \n');  
        fprintf('                                                            \n');
        fprintf('A          = the input matrix                               \n');
        fprintf('                                                            \n');
        fprintf('m          = subspace dimension                             \n');
        fprintf('                                                            \n');
        fprintf('ver        = method                                         \n');
        fprintf('                                                            \n');
        fprintf('percentage = the percentage of the trace to be retained     \n');
        fprintf('                                                            \n');
        fprintf('eps        = threshold value to stop the subspace iteration \n');
        fprintf('                                                            \n');
        fprintf('maxit      = the maximum number of iterations               \n');
        fprintf('                                                            \n');
        fprintf('                                                            \n');
        fprintf('                                                            \n');
        fprintf('OUTPUT arguments:                                           \n');
        fprintf('                                                            \n');
        fprintf('V          = a matrix containing eigenvector                \n');
        fprintf('                                                            \n');
        fprintf('w          = a vector containing the computed eigenvalues   \n');
        fprintf('                                                            \n');
        fprintf('res_ev     = a vector containing the norm of the residual   \n');
        fprintf('             for each eigenvalue                            \n');
        fprintf('                                                            \n');
        fprintf('it_ev      = a vector containing the number of iterations   \n');
        fprintf('             for each eigenvalue                            \n');
        fprintf(' \n');
        V=0;
        w=0;
	res_ev = 0;
	it_ev = 0;
        return
    else
        [V, w, res_ev, it_ev]  = mex_subspace_iter_ev(A, m, ver, percentage, eps, maxit);
    end
