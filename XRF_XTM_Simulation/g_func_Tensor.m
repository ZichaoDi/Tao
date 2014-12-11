% Generated by ADiMat 0.6.0-4728
% Copyright 2009-2013 Johannes Willkomm, Fachgebiet Scientific Computing,
% TU Darmstadt, 64289 Darmstadt, Germany
% Copyright 2001-2008 Andre Vehreschild, Institute for Scientific Computing,
% RWTH Aachen University, 52056 Aachen, Germany.
% Visit us on the web at http://www.adimat.de
% Report bugs to adimat-users@lists.sc.informatik.tu-darmstadt.de
%
%
%                             DISCLAIMER
%
% ADiMat was prepared as part of an employment at the Institute
% for Scientific Computing, RWTH Aachen University, Germany and is
% provided AS IS. NEITHER THE AUTHOR(S), THE GOVERNMENT OF THE FEDERAL
% REPUBLIC OF GERMANY NOR ANY AGENCY THEREOF, NOR THE RWTH AACHEN UNIVERSITY,
% INCLUDING ANY OF THEIR EMPLOYEES OR OFFICERS, MAKES ANY WARRANTY,
% EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY
% FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION OR
% PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
% PRIVATELY OWNED RIGHTS.
%
% Global flags were:
% FORWARDMODE -- Apply the forward mode to the files.
% NOOPEROPTIM -- Do not use optimized operators. I.e.:
%		 g_a*b*g_c -/-> mtimes3(g_a, b, g_c)
% NOLOCALCSE  -- Do not use local common subexpression elimination when
%		 canonicalizing the code.
% NOGLOBALCSE -- Prevents the application of global common subexpression
%		 elimination after canonicalizing the code.
% NOPRESCALARFOLDING -- Switch off folding of scalar constants before
%		 augmentation.
% NOPOSTSCALARFOLDING -- Switch off folding of scalar constants after
%		 augmentation.
% NOCONSTFOLDMULT0 -- Switch off folding of product with one factor
%		 being zero: b*0=0.
% FUNCMODE    -- Inputfile is a function (This flag can not be set explicitly).
% NOTMPCLEAR  -- Suppress generation of clear g_* instructions.
% UNBOUND_XML  -- Write list of unbound identifiers in XML format.
% DEPENDENCIES_XML  -- Write list of functions in XML format.
% UNBOUND_ERROR	-- Stop with error if unbound identifiers found (default).
% FUNCTION_LIST_XML	-- Write list of functions to XML file.
% VERBOSITYLEVEL=5
% AD_IVARS= W
% AD_DVARS= f

function [g_f, f]= g_func_Tensor(g_W, W, xrfData, M, NumElement, L, GlobalInd, SelfInd, thetan, m, nTau)
   global NumSSDlet numChannel NoSelfAbsorption
   global mtol
   f= 0; 
   g_f= g_zeros(size(f));
   g_tmp_func_Tensor_00011= call(@reshape, g_W, mtol, NumElement);
   tmp_func_Tensor_00011= reshape(W, mtol, NumElement); 
   % Update detected: W= some_expression(W,...)
   g_W= g_tmp_func_Tensor_00011;
   W= tmp_func_Tensor_00011;
   tmp_func_Tensor_00012= reshape(L, length(thetan), nTau+ 1, mtol); 
   %%%%% ====================================================================
   % Update detected: L= some_expression(L,...)
   L= tmp_func_Tensor_00012;
   for n= 1: length(thetan)
      sum_Tau= 0; 
      g_sum_Tau= g_zeros(size(sum_Tau));
      tmp_func_Tensor_00000= nTau+ 1;
      for i= 1: tmp_func_Tensor_00000
         XRF_v= zeros(1, numChannel); 
         g_XRF_v= g_zeros(size(XRF_v));
         index= GlobalInd{n, i}; 
         if (~isempty(index))
            index_sub= sub2ind(m, index(: , 2), index(: , 1)); 
            msub= length(index_sub); 
            for v_count= 1: msub
               v= index_sub(v_count); 
               if (isempty(SelfInd{n, i, v}{1}))
                  InTens= 1; 
                  g_InTens= g_zeros(size(InTens));
               else 
                  g_tmp_W_00000= g_W(SelfInd{n, i, v}{1}, : );
                  tmp_W_00000= W(SelfInd{n, i, v}{1}, : );
                  g_tmp_func_Tensor_00001= g_tmp_W_00000.* SelfInd{n, i, v}{3};
                  tmp_func_Tensor_00001= tmp_W_00000.* SelfInd{n, i, v}{3};
                  g_tmp_sum_00000= call(@sum, g_tmp_func_Tensor_00001);
                  tmp_sum_00000= sum(tmp_func_Tensor_00001);
                  g_tmp_sum_00001= call(@sum, g_tmp_sum_00000);
                  tmp_sum_00001= sum(tmp_sum_00000);
                  g_tmp_func_Tensor_00002= -g_tmp_sum_00001;
                  tmp_func_Tensor_00002= -tmp_sum_00001;
                  g_InTens= g_tmp_func_Tensor_00002.* exp(tmp_func_Tensor_00002);
                  InTens= exp(tmp_func_Tensor_00002); 
               end
               
               if (isempty(SelfInd{n, i, v}{2})| NoSelfAbsorption)
                  OutTens_d= ones(1, NumElement); 
                  g_OutTens_d= g_zeros(size(OutTens_d));
               else 
                  OutTens_d= 0; 
                  g_OutTens_d= g_zeros(size(OutTens_d));
                  for d= 1: NumSSDlet%                         OutTens_d=OutTens_d+exp(-sum(sum(bsxfun(@times,W(SelfInd{n,i,v}{2}{d},:),reshape(SelfInd{n,i,v}{4}{d},length(SelfInd{n,i,v}{2}{d}),NumElement,NumElement)),1),2));
                     g_tmp_W_00001= g_W(SelfInd{n, i, v}{2}{d}, : );
                     tmp_W_00001= W(SelfInd{n, i, v}{2}{d}, : );
                     tmp_func_Tensor_00003= [1, 1, NumElement];
                     g_tmp_repmat_00000= call(@repmat, g_tmp_W_00001, tmp_func_Tensor_00003);
                     tmp_repmat_00000= repmat(tmp_W_00001, tmp_func_Tensor_00003);
                     g_Tmmp= g_tmp_repmat_00000.* reshape(SelfInd{n, i, v}{4}{d}, length(SelfInd{n, i, v}{2}{d}), NumElement, NumElement);
                     Tmmp= tmp_repmat_00000.* reshape(SelfInd{n, i, v}{4}{d}, length(SelfInd{n, i, v}{2}{d}), NumElement, NumElement); 
                     g_tmp_sum_00002= call(@sum, g_Tmmp, 1);
                     tmp_sum_00002= sum(Tmmp, 1);
                     g_tmp_sum_00003= call(@sum, g_tmp_sum_00002, 2);
                     tmp_sum_00003= sum(tmp_sum_00002, 2);
                     g_tmp_func_Tensor_00004= -g_tmp_sum_00003;
                     tmp_func_Tensor_00004= -tmp_sum_00003;
                     g_tmp_exp_00000= g_tmp_func_Tensor_00004.* exp(tmp_func_Tensor_00004);
                     tmp_exp_00000= exp(tmp_func_Tensor_00004);
                     g_tmp_func_Tensor_00013= g_OutTens_d+ g_tmp_exp_00000;
                     tmp_func_Tensor_00013= OutTens_d+ tmp_exp_00000; 
                     % Update detected: OutTens_d= some_expression(OutTens_d,...)
                     g_OutTens_d= g_tmp_func_Tensor_00013;
                     OutTens_d= tmp_func_Tensor_00013;
                  end
                  g_tmp_func_Tensor_00014= adimat_g_mrdivide1(g_OutTens_d, OutTens_d, NumSSDlet);
                  tmp_func_Tensor_00014= OutTens_d/ NumSSDlet; 
                  % Update detected: OutTens_d= some_expression(OutTens_d,...)
                  g_OutTens_d= g_tmp_func_Tensor_00014;
                  OutTens_d= tmp_func_Tensor_00014;
               end
               g_tmp_reshape_00000= call(@reshape, g_OutTens_d, 1, NumElement);
               tmp_reshape_00000= reshape(OutTens_d, 1, NumElement);
               g_tmp_func_Tensor_00005= g_InTens* tmp_reshape_00000+ InTens* g_tmp_reshape_00000;
               tmp_func_Tensor_00005= InTens* tmp_reshape_00000;
               g_tmp_W_00002= g_W(v, : );
               tmp_W_00002= W(v, : );
               g_tmp_func_Tensor_00006= g_tmp_func_Tensor_00005.* tmp_W_00002+ tmp_func_Tensor_00005.* g_tmp_W_00002;
               tmp_func_Tensor_00006= tmp_func_Tensor_00005.* tmp_W_00002;
               g_tmp_func_Tensor_00007= L(n, i, v)* g_tmp_func_Tensor_00006* M;
               tmp_func_Tensor_00007= L(n, i, v)* tmp_func_Tensor_00006* M;
               g_tmp_func_Tensor_00015= g_XRF_v+ g_tmp_func_Tensor_00007;
               tmp_func_Tensor_00015= XRF_v+ tmp_func_Tensor_00007; 
               % Update detected: XRF_v= some_expression(XRF_v,...)
               g_XRF_v= g_tmp_func_Tensor_00015;
               XRF_v= tmp_func_Tensor_00015;
            end
            g_tmp_func_Tensor_00008= -g_XRF_v+ g_zeros(size(xrfData{n, i}));
            tmp_func_Tensor_00008= xrfData{n, i}- XRF_v;
            g_tmp_func_Tensor_00009= -g_XRF_v+ g_zeros(size(xrfData{n, i}));
            tmp_func_Tensor_00009= xrfData{n, i}- XRF_v;
            g_tmp_func_Tensor_00010= g_tmp_func_Tensor_00008* (tmp_func_Tensor_00009)' + tmp_func_Tensor_00008* (g_tmp_func_Tensor_00009)' ;
            tmp_func_Tensor_00010= tmp_func_Tensor_00008* (tmp_func_Tensor_00009)' ;
            g_tmp_func_Tensor_00016= g_sum_Tau+ g_tmp_func_Tensor_00010;
            tmp_func_Tensor_00016= sum_Tau+ tmp_func_Tensor_00010; 
            % Update detected: sum_Tau= some_expression(sum_Tau,...)
            g_sum_Tau= g_tmp_func_Tensor_00016;
            sum_Tau= tmp_func_Tensor_00016;
         end
         tmp_func_Tensor_00000= nTau+ 1;
      end
      g_tmp_func_Tensor_00017= g_f+ g_sum_Tau;
      tmp_func_Tensor_00017= f+ sum_Tau; 
      % Update detected: f= some_expression(f,...)
      g_f= g_tmp_func_Tensor_00017;
      f= tmp_func_Tensor_00017;
   end
end