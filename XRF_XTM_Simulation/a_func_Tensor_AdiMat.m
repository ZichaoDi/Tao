% Generated by ADiMat 0.6.0-4728
% © 2001-2008 Andre Vehreschild <vehreschild@sc.rwth-aachen.de>
% © 2009-2013 Johannes Willkomm <johannes.willkomm@sc.tu-darmstadt.de>
% RWTH Aachen University, 52056 Aachen, Germany
% TU Darmstadt, 64289 Darmstadt, Germany
% Visit us on the web at http://www.adimat.de/
% Report bugs to adimat-users@lists.sc.informatik.tu-darmstadt.de
%
%                             DISCLAIMER
% 
% ADiMat was prepared as part of an employment at the Institute for Scientific Computing,
% RWTH Aachen University, Germany and at the Institute for Scientific Computing,
% TU Darmstadt, Germany and is provided AS IS. 
% NEITHER THE AUTHOR(S), THE GOVERNMENT OF THE FEDERAL REPUBLIC OF GERMANY
% NOR ANY AGENCY THEREOF, NOR THE RWTH AACHEN UNIVERSITY, NOT THE TU DARMSTADT,
% INCLUDING ANY OF THEIR EMPLOYEES OR OFFICERS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
% OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,
% OR USEFULNESS OF ANY INFORMATION OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE
% WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
%
% Flags: BACKWARDMODE,  NOOPEROPTIM,
%   NOLOCALCSE,  NOGLOBALCSE,  NOPRESCALARFOLDING,
%   NOPOSTSCALARFOLDING,  NOCONSTFOLDMULT0,  FUNCMODE,
%   NOTMPCLEAR,  DUMP_XML,  UNBOUND_XML,
%   DEPENDENCIES_XML,  PARSE_ONLY,  UNBOUND_ERROR,
%   FUNCTION_LIST_XML
%
% Parameters:
%  - dependents=f
%  - independents=W
%  - inputEncoding=ISO-8859-1
%  - output-mode: plain
%  - output-file: outputs/a_func_Tensor_AdiMat.m
%  - output-file-prefix: 
%  - output-directory: outputs
%
% Functions in this file: a_func_Tensor_AdiMat, rec_func_Tensor_AdiMat,
%  ret_func_Tensor_AdiMat
%

function [a_W nr_f] = a_func_Tensor_AdiMat(W, xrfData, M, NumElement, L, GlobalInd, SelfInd, thetan, m, nTau, a_f)
   global NumSSDlet numChannel NoSelfAbsorption mtol;
   tmpca3 = 0;
   tmpca2 = 0;
   tmpca1 = 0;
   tmpda4 = 0;
   tmpda3 = 0;
   tmpda2 = 0;
   tmpca5 = 0;
   tmpca4 = 0;
   sum_Tau = 0;
   XRF_v = 0;
   index = 0;
   index_sub = 0;
   msub = 0;
   v = 0;
   InTens = 0;
   OutTens_d = 0;
   Tmmp = 0;
   f = 0;
   adimat_push1(W);
   W = reshape(W, mtol, NumElement);
   adimat_push1(L);
   L = reshape(L, length(thetan), nTau + 1, mtol);
   tmpfra1_2 = length(thetan);
   for n=1 : tmpfra1_2
      adimat_push1(sum_Tau);
      sum_Tau = 0;
      tmpfra2_2 = nTau + 1;
      for i=1 : tmpfra2_2
         adimat_push1(XRF_v);
         XRF_v = zeros(1, numChannel);
         adimat_push1(index);
         index = GlobalInd{n, i};
         tmpba1 = 0;
         if ~isempty(index)
            tmpba1 = 1;
            adimat_push1(index_sub);
            index_sub = sub2ind(m, index(:, 2), index(:, 1));
            adimat_push1(msub);
            msub = length(index_sub);
            tmpfra3_2 = msub;
            for v_count=1 : tmpfra3_2
               adimat_push1(v);
               v = index_sub(v_count);
               tmpba2 = 0;
               if isempty(SelfInd{n, i, v}{1})
                  tmpba2 = 1;
                  adimat_push1(InTens);
                  InTens = 1;
               else
                  adimat_push1(tmpca3);
                  tmpca3 = W(SelfInd{n, i, v}{1}, :) .* SelfInd{n, i, v}{3};
                  adimat_push1(tmpca2);
                  tmpca2 = sum(tmpca3);
                  adimat_push1(tmpca1);
                  tmpca1 = sum(tmpca2);
                  adimat_push1(InTens);
                  InTens = exp(-tmpca1);
               end
               adimat_push1(tmpba2);
               tmpba2 = 0;
               if isempty(SelfInd{n, i, v}{2}) | NoSelfAbsorption
                  tmpba2 = 1;
                  adimat_push1(OutTens_d);
                  OutTens_d = ones(1, NumElement);
               else
                  adimat_push1(OutTens_d);
                  OutTens_d = zeros(1, 1, NumElement);
                  tmpfra4_2 = NumSSDlet;
                  for d=1 : tmpfra4_2
                     tmpba3 = 0;
                     if ~isempty(SelfInd{n, i, v}{2}{d})
                        tmpba3 = 1;
                        adimat_push1(tmpda4);
                        tmpda4 = length(SelfInd{n, i, v}{2}{d});
                        adimat_push1(tmpda3);
                        tmpda3 = reshape(SelfInd{n, i, v}{4}{d}, tmpda4, NumElement, NumElement);
                        adimat_push1(tmpda2);
                        tmpda2 = [1 1 NumElement];
                        adimat_push1(tmpca1);
                        tmpca1 = repmat(W(SelfInd{n, i, v}{2}{d}, :), tmpda2);
                        adimat_push1(Tmmp);
                        Tmmp = tmpca1 .* tmpda3;
                     else
                        adimat_push1(Tmmp);
                        Tmmp = 0;
                     end
                     adimat_push(tmpba3, tmpca3);
                     tmpca3 = sum(Tmmp, 1);
                     adimat_push1(tmpca2);
                     tmpca2 = sum(tmpca3, 2);
                     adimat_push1(tmpca1);
                     tmpca1 = exp(-tmpca2);
                     adimat_push1(OutTens_d);
                     OutTens_d = OutTens_d + tmpca1;
                  end
                  adimat_push(tmpfra4_2, OutTens_d);
                  OutTens_d = OutTens_d / NumSSDlet;
               end
               adimat_push(tmpba2, tmpca5);
               tmpca5 = reshape(OutTens_d, 1, NumElement);
               adimat_push1(tmpca4);
               tmpca4 = InTens * tmpca5;
               adimat_push1(tmpca3);
               tmpca3 = tmpca4 .* W(v, :);
               adimat_push1(tmpca2);
               tmpca2 = L(n, i, v) * tmpca3;
               adimat_push1(tmpca1);
               tmpca1 = tmpca2 * M;
               adimat_push1(XRF_v);
               XRF_v = XRF_v + tmpca1;
            end
            adimat_push(tmpfra3_2, tmpca3);
            tmpca3 = xrfData{n, i} - XRF_v;
            adimat_push1(tmpca2);
            tmpca2 = xrfData{n, i} - XRF_v;
            adimat_push1(tmpca1);
            tmpca1 = tmpca2 * tmpca3';
            adimat_push1(sum_Tau);
            sum_Tau = sum_Tau + tmpca1;
         end
         adimat_push1(tmpba1);
      end
      adimat_push(tmpfra2_2, f);
      f = f + sum_Tau;
   end
   adimat_push1(tmpfra1_2);
   nr_f = f;
   [a_sum_Tau a_XRF_v a_InTens a_OutTens_d a_Tmmp a_tmpca3 a_tmpca2 a_tmpca1 a_tmpca5 a_tmpca4 a_W] = a_zeros(sum_Tau, XRF_v, InTens, OutTens_d, Tmmp, tmpca3, tmpca2, tmpca1, tmpca5, tmpca4, W);
   if nargin < 11
      a_f = a_zeros1(f);
   end
   tmpfra1_2 = adimat_pop1;
   for n=fliplr(1 : tmpfra1_2)
      f = adimat_pop1;
      a_sum_Tau = adimat_adjsum(a_sum_Tau, adimat_adjred(sum_Tau, a_f));
      tmpsa1 = a_f;
      a_f = a_zeros1(f);
      a_f = adimat_adjsum(a_f, adimat_adjred(f, tmpsa1));
      tmpfra2_2 = adimat_pop1;
      for i=fliplr(1 : tmpfra2_2)
         tmpba1 = adimat_pop1;
         if tmpba1 == 1
            sum_Tau = adimat_pop1;
            a_tmpca1 = adimat_adjsum(a_tmpca1, adimat_adjred(tmpca1, a_sum_Tau));
            tmpsa1 = a_sum_Tau;
            a_sum_Tau = a_zeros1(sum_Tau);
            a_sum_Tau = adimat_adjsum(a_sum_Tau, adimat_adjred(sum_Tau, tmpsa1));
            tmpca1 = adimat_pop1;
            a_tmpca2 = adimat_adjsum(a_tmpca2, adimat_adjmultl(tmpca2, a_tmpca1, tmpca3'));
            a_tmpca3 = adimat_adjsum(a_tmpca3, a_ctranspose(adimat_adjmultr(tmpca3', tmpca2, a_tmpca1), tmpca3));
            a_tmpca1 = a_zeros1(tmpca1);
            tmpca2 = adimat_pop1;
            a_XRF_v = adimat_adjsum(a_XRF_v, adimat_adjred(XRF_v, -a_tmpca2));
            a_tmpca2 = a_zeros1(tmpca2);
            tmpca3 = adimat_pop1;
            a_XRF_v = adimat_adjsum(a_XRF_v, adimat_adjred(XRF_v, -a_tmpca3));
            a_tmpca3 = a_zeros1(tmpca3);
            tmpfra3_2 = adimat_pop1;
            for v_count=fliplr(1 : tmpfra3_2)
               XRF_v = adimat_pop1;
               a_tmpca1 = adimat_adjsum(a_tmpca1, adimat_adjred(tmpca1, a_XRF_v));
               tmpsa1 = a_XRF_v;
               a_XRF_v = a_zeros1(XRF_v);
               a_XRF_v = adimat_adjsum(a_XRF_v, adimat_adjred(XRF_v, tmpsa1));
               tmpca1 = adimat_pop1;
               a_tmpca2 = adimat_adjsum(a_tmpca2, adimat_adjmultl(tmpca2, a_tmpca1, M));
               a_tmpca1 = a_zeros1(tmpca1);
               tmpca2 = adimat_pop1;
               a_tmpca3 = adimat_adjsum(a_tmpca3, adimat_adjmultr(tmpca3, L(n, i, v), a_tmpca2));
               a_tmpca2 = a_zeros1(tmpca2);
               tmpca3 = adimat_pop1;
               a_tmpca4 = adimat_adjsum(a_tmpca4, adimat_adjred(tmpca4, a_tmpca3 .* W(v, :)));
               a_W(v, :) = adimat_adjsum(a_W(v, :), adimat_adjred(W(v, :), tmpca4 .* a_tmpca3));
               a_tmpca3 = a_zeros1(tmpca3);
               tmpca4 = adimat_pop1;
               a_InTens = adimat_adjsum(a_InTens, adimat_adjmultl(InTens, a_tmpca4, tmpca5));
               a_tmpca5 = adimat_adjsum(a_tmpca5, adimat_adjmultr(tmpca5, InTens, a_tmpca4));
               a_tmpca4 = a_zeros1(tmpca4);
               tmpca5 = adimat_pop1;
               a_OutTens_d = adimat_adjsum(a_OutTens_d, reshape(a_tmpca5, size(OutTens_d)));
               a_tmpca5 = a_zeros1(tmpca5);
               tmpba2 = adimat_pop1;
               if tmpba2 == 1
                  OutTens_d = adimat_pop1;
                  a_OutTens_d = a_zeros1(OutTens_d);
               else
                  [tmpadjc1] = adimat_a_mrdividel(OutTens_d, NumSSDlet, a_OutTens_d);
                  OutTens_d = adimat_pop1;
                  tmpsa1 = a_OutTens_d;
                  a_OutTens_d = a_zeros1(OutTens_d);
                  a_OutTens_d = adimat_adjsum(a_OutTens_d, tmpadjc1);
                  tmpfra4_2 = adimat_pop1;
                  for d=fliplr(1 : tmpfra4_2)
                     OutTens_d = adimat_pop1;
                     a_tmpca1 = adimat_adjsum(a_tmpca1, adimat_adjred(tmpca1, a_OutTens_d));
                     tmpsa1 = a_OutTens_d;
                     a_OutTens_d = a_zeros1(OutTens_d);
                     a_OutTens_d = adimat_adjsum(a_OutTens_d, adimat_adjred(OutTens_d, tmpsa1));
                     tmpca1 = adimat_pop1;
                     a_tmpca2 = adimat_adjsum(a_tmpca2, -(exp(-tmpca2) .* a_tmpca1));
                     a_tmpca1 = a_zeros1(tmpca1);
                     tmpca2 = adimat_pop1;
                     a_tmpca3 = adimat_adjsum(a_tmpca3, a_sum(a_tmpca2, tmpca3, 2));
                     a_tmpca2 = a_zeros1(tmpca2);
                     tmpca3 = adimat_pop1;
                     a_Tmmp = adimat_adjsum(a_Tmmp, a_sum(a_tmpca3, Tmmp, 1));
                     a_tmpca3 = a_zeros1(tmpca3);
                     tmpba3 = adimat_pop1;
                     if tmpba3 == 1
                        Tmmp = adimat_pop1;
                        a_tmpca1 = adimat_adjsum(a_tmpca1, adimat_adjred(tmpca1, a_Tmmp .* tmpda3));
                        a_Tmmp = a_zeros1(Tmmp);
                        tmpca1 = adimat_pop1;
                        a_W(SelfInd{n, i, v}{2}{d}, :) = adimat_adjsum(a_W(SelfInd{n, i, v}{2}{d}, :), a_repmat(a_tmpca1, W(SelfInd{n, i, v}{2}{d}, :), tmpda2));
                        a_tmpca1 = a_zeros1(tmpca1);
                        [tmpda2 tmpda3 tmpda4] = adimat_pop;
                     else
                        Tmmp = adimat_pop1;
                        a_Tmmp = a_zeros1(Tmmp);
                     end
                  end
                  OutTens_d = adimat_pop1;
                  a_OutTens_d = a_zeros1(OutTens_d);
               end
               tmpba2 = adimat_pop1;
               if tmpba2 == 1
                  InTens = adimat_pop1;
                  a_InTens = a_zeros1(InTens);
               else
                  InTens = adimat_pop1;
                  a_tmpca1 = adimat_adjsum(a_tmpca1, -(exp(-tmpca1) .* a_InTens));
                  a_InTens = a_zeros1(InTens);
                  tmpca1 = adimat_pop1;
                  a_tmpca2 = adimat_adjsum(a_tmpca2, a_sum(a_tmpca1, tmpca2));
                  a_tmpca1 = a_zeros1(tmpca1);
                  tmpca2 = adimat_pop1;
                  a_tmpca3 = adimat_adjsum(a_tmpca3, a_sum(a_tmpca2, tmpca3));
                  a_tmpca2 = a_zeros1(tmpca2);
                  tmpca3 = adimat_pop1;
                  a_W(SelfInd{n, i, v}{1}, :) = adimat_adjsum(a_W(SelfInd{n, i, v}{1}, :), adimat_adjred(W(SelfInd{n, i, v}{1}, :), a_tmpca3 .* SelfInd{n, i, v}{3}));
                  a_tmpca3 = a_zeros1(tmpca3);
               end
               v = adimat_pop1;
            end
            [msub index_sub] = adimat_pop;
         end
         [index XRF_v] = adimat_pop;
         a_XRF_v = a_zeros1(XRF_v);
      end
      sum_Tau = adimat_pop1;
      a_sum_Tau = a_zeros1(sum_Tau);
   end
   [L W] = adimat_pop;
   tmpsa1 = a_W;
   a_W = a_zeros1(W);
   a_W = adimat_adjsum(a_W, reshape(tmpsa1, size(W)));
end

function f = rec_func_Tensor_AdiMat(W, xrfData, M, NumElement, L, GlobalInd, SelfInd, thetan, m, nTau)
   global NumSSDlet numChannel NoSelfAbsorption mtol;
   tmpca3 = 0;
   tmpca2 = 0;
   tmpca1 = 0;
   tmpda4 = 0;
   tmpda3 = 0;
   tmpda2 = 0;
   tmpca5 = 0;
   tmpca4 = 0;
   sum_Tau = 0;
   XRF_v = 0;
   index = 0;
   index_sub = 0;
   msub = 0;
   v = 0;
   InTens = 0;
   OutTens_d = 0;
   Tmmp = 0;
   f = 0;
   adimat_push1(W);
   W = reshape(W, mtol, NumElement);
   adimat_push1(L);
   L = reshape(L, length(thetan), nTau + 1, mtol);
   tmpfra1_2 = length(thetan);
   for n=1 : tmpfra1_2
      adimat_push1(sum_Tau);
      sum_Tau = 0;
      tmpfra2_2 = nTau + 1;
      for i=1 : tmpfra2_2
         adimat_push1(XRF_v);
         XRF_v = zeros(1, numChannel);
         adimat_push1(index);
         index = GlobalInd{n, i};
         tmpba1 = 0;
         if ~isempty(index)
            tmpba1 = 1;
            adimat_push1(index_sub);
            index_sub = sub2ind(m, index(:, 2), index(:, 1));
            adimat_push1(msub);
            msub = length(index_sub);
            tmpfra3_2 = msub;
            for v_count=1 : tmpfra3_2
               adimat_push1(v);
               v = index_sub(v_count);
               tmpba2 = 0;
               if isempty(SelfInd{n, i, v}{1})
                  tmpba2 = 1;
                  adimat_push1(InTens);
                  InTens = 1;
               else
                  adimat_push1(tmpca3);
                  tmpca3 = W(SelfInd{n, i, v}{1}, :) .* SelfInd{n, i, v}{3};
                  adimat_push1(tmpca2);
                  tmpca2 = sum(tmpca3);
                  adimat_push1(tmpca1);
                  tmpca1 = sum(tmpca2);
                  adimat_push1(InTens);
                  InTens = exp(-tmpca1);
               end
               adimat_push1(tmpba2);
               tmpba2 = 0;
               if isempty(SelfInd{n, i, v}{2}) | NoSelfAbsorption
                  tmpba2 = 1;
                  adimat_push1(OutTens_d);
                  OutTens_d = ones(1, NumElement);
               else
                  adimat_push1(OutTens_d);
                  OutTens_d = zeros(1, 1, NumElement);
                  tmpfra4_2 = NumSSDlet;
                  for d=1 : tmpfra4_2
                     tmpba3 = 0;
                     if ~isempty(SelfInd{n, i, v}{2}{d})
                        tmpba3 = 1;
                        adimat_push1(tmpda4);
                        tmpda4 = length(SelfInd{n, i, v}{2}{d});
                        adimat_push1(tmpda3);
                        tmpda3 = reshape(SelfInd{n, i, v}{4}{d}, tmpda4, NumElement, NumElement);
                        adimat_push1(tmpda2);
                        tmpda2 = [1 1 NumElement];
                        adimat_push1(tmpca1);
                        tmpca1 = repmat(W(SelfInd{n, i, v}{2}{d}, :), tmpda2);
                        adimat_push1(Tmmp);
                        Tmmp = tmpca1 .* tmpda3;
                     else
                        adimat_push1(Tmmp);
                        Tmmp = 0;
                     end
                     adimat_push(tmpba3, tmpca3);
                     tmpca3 = sum(Tmmp, 1);
                     adimat_push1(tmpca2);
                     tmpca2 = sum(tmpca3, 2);
                     adimat_push1(tmpca1);
                     tmpca1 = exp(-tmpca2);
                     adimat_push1(OutTens_d);
                     OutTens_d = OutTens_d + tmpca1;
                  end
                  adimat_push(tmpfra4_2, OutTens_d);
                  OutTens_d = OutTens_d / NumSSDlet;
               end
               adimat_push(tmpba2, tmpca5);
               tmpca5 = reshape(OutTens_d, 1, NumElement);
               adimat_push1(tmpca4);
               tmpca4 = InTens * tmpca5;
               adimat_push1(tmpca3);
               tmpca3 = tmpca4 .* W(v, :);
               adimat_push1(tmpca2);
               tmpca2 = L(n, i, v) * tmpca3;
               adimat_push1(tmpca1);
               tmpca1 = tmpca2 * M;
               adimat_push1(XRF_v);
               XRF_v = XRF_v + tmpca1;
            end
            adimat_push(tmpfra3_2, tmpca3);
            tmpca3 = xrfData{n, i} - XRF_v;
            adimat_push1(tmpca2);
            tmpca2 = xrfData{n, i} - XRF_v;
            adimat_push1(tmpca1);
            tmpca1 = tmpca2 * tmpca3';
            adimat_push1(sum_Tau);
            sum_Tau = sum_Tau + tmpca1;
         end
         adimat_push1(tmpba1);
      end
      adimat_push(tmpfra2_2, f);
      f = f + sum_Tau;
   end
   adimat_push(tmpfra1_2, sum_Tau, XRF_v, index, index_sub, msub, v, InTens, OutTens_d, Tmmp, tmpca3, tmpca2, tmpca1, tmpda4, tmpda3, tmpda2, tmpca5, tmpca4, f, W);
   if nargin > 1
      adimat_push1(xrfData);
   end
   if nargin > 2
      adimat_push1(M);
   end
   if nargin > 3
      adimat_push1(NumElement);
   end
   if nargin > 4
      adimat_push1(L);
   end
   if nargin > 5
      adimat_push1(GlobalInd);
   end
   if nargin > 6
      adimat_push1(SelfInd);
   end
   if nargin > 7
      adimat_push1(thetan);
   end
   if nargin > 8
      adimat_push1(m);
   end
   if nargin > 9
      adimat_push1(nTau);
   end
   adimat_push1(nargin);
end

function a_W = ret_func_Tensor_AdiMat(a_f)
   global NumSSDlet numChannel NoSelfAbsorption mtol;
   tmpnargin = adimat_pop1;
   if tmpnargin > 9
      nTau = adimat_pop1;
   end
   if tmpnargin > 8
      m = adimat_pop1;
   end
   if tmpnargin > 7
      thetan = adimat_pop1;
   end
   if tmpnargin > 6
      SelfInd = adimat_pop1;
   end
   if tmpnargin > 5
      GlobalInd = adimat_pop1;
   end
   if tmpnargin > 4
      L = adimat_pop1;
   end
   if tmpnargin > 3
      NumElement = adimat_pop1;
   end
   if tmpnargin > 2
      M = adimat_pop1;
   end
   if tmpnargin > 1
      xrfData = adimat_pop1;
   end
   [W f tmpca4 tmpca5 tmpda2 tmpda3 tmpda4 tmpca1 tmpca2 tmpca3 Tmmp OutTens_d InTens v msub index_sub index XRF_v sum_Tau] = adimat_pop;
   [a_sum_Tau a_XRF_v a_InTens a_OutTens_d a_Tmmp a_tmpca3 a_tmpca2 a_tmpca1 a_tmpca5 a_tmpca4 a_W] = a_zeros(sum_Tau, XRF_v, InTens, OutTens_d, Tmmp, tmpca3, tmpca2, tmpca1, tmpca5, tmpca4, W);
   if nargin < 1
      a_f = a_zeros1(f);
   end
   tmpfra1_2 = adimat_pop1;
   for n=fliplr(1 : tmpfra1_2)
      f = adimat_pop1;
      a_sum_Tau = adimat_adjsum(a_sum_Tau, adimat_adjred(sum_Tau, a_f));
      tmpsa1 = a_f;
      a_f = a_zeros1(f);
      a_f = adimat_adjsum(a_f, adimat_adjred(f, tmpsa1));
      tmpfra2_2 = adimat_pop1;
      for i=fliplr(1 : tmpfra2_2)
         tmpba1 = adimat_pop1;
         if tmpba1 == 1
            sum_Tau = adimat_pop1;
            a_tmpca1 = adimat_adjsum(a_tmpca1, adimat_adjred(tmpca1, a_sum_Tau));
            tmpsa1 = a_sum_Tau;
            a_sum_Tau = a_zeros1(sum_Tau);
            a_sum_Tau = adimat_adjsum(a_sum_Tau, adimat_adjred(sum_Tau, tmpsa1));
            tmpca1 = adimat_pop1;
            a_tmpca2 = adimat_adjsum(a_tmpca2, adimat_adjmultl(tmpca2, a_tmpca1, tmpca3'));
            a_tmpca3 = adimat_adjsum(a_tmpca3, a_ctranspose(adimat_adjmultr(tmpca3', tmpca2, a_tmpca1), tmpca3));
            a_tmpca1 = a_zeros1(tmpca1);
            tmpca2 = adimat_pop1;
            a_XRF_v = adimat_adjsum(a_XRF_v, adimat_adjred(XRF_v, -a_tmpca2));
            a_tmpca2 = a_zeros1(tmpca2);
            tmpca3 = adimat_pop1;
            a_XRF_v = adimat_adjsum(a_XRF_v, adimat_adjred(XRF_v, -a_tmpca3));
            a_tmpca3 = a_zeros1(tmpca3);
            tmpfra3_2 = adimat_pop1;
            for v_count=fliplr(1 : tmpfra3_2)
               XRF_v = adimat_pop1;
               a_tmpca1 = adimat_adjsum(a_tmpca1, adimat_adjred(tmpca1, a_XRF_v));
               tmpsa1 = a_XRF_v;
               a_XRF_v = a_zeros1(XRF_v);
               a_XRF_v = adimat_adjsum(a_XRF_v, adimat_adjred(XRF_v, tmpsa1));
               tmpca1 = adimat_pop1;
               a_tmpca2 = adimat_adjsum(a_tmpca2, adimat_adjmultl(tmpca2, a_tmpca1, M));
               a_tmpca1 = a_zeros1(tmpca1);
               tmpca2 = adimat_pop1;
               a_tmpca3 = adimat_adjsum(a_tmpca3, adimat_adjmultr(tmpca3, L(n, i, v), a_tmpca2));
               a_tmpca2 = a_zeros1(tmpca2);
               tmpca3 = adimat_pop1;
               a_tmpca4 = adimat_adjsum(a_tmpca4, adimat_adjred(tmpca4, a_tmpca3 .* W(v, :)));
               a_W(v, :) = adimat_adjsum(a_W(v, :), adimat_adjred(W(v, :), tmpca4 .* a_tmpca3));
               a_tmpca3 = a_zeros1(tmpca3);
               tmpca4 = adimat_pop1;
               a_InTens = adimat_adjsum(a_InTens, adimat_adjmultl(InTens, a_tmpca4, tmpca5));
               a_tmpca5 = adimat_adjsum(a_tmpca5, adimat_adjmultr(tmpca5, InTens, a_tmpca4));
               a_tmpca4 = a_zeros1(tmpca4);
               tmpca5 = adimat_pop1;
               a_OutTens_d = adimat_adjsum(a_OutTens_d, reshape(a_tmpca5, size(OutTens_d)));
               a_tmpca5 = a_zeros1(tmpca5);
               tmpba2 = adimat_pop1;
               if tmpba2 == 1
                  OutTens_d = adimat_pop1;
                  a_OutTens_d = a_zeros1(OutTens_d);
               else
                  [tmpadjc1] = adimat_a_mrdividel(OutTens_d, NumSSDlet, a_OutTens_d);
                  OutTens_d = adimat_pop1;
                  tmpsa1 = a_OutTens_d;
                  a_OutTens_d = a_zeros1(OutTens_d);
                  a_OutTens_d = adimat_adjsum(a_OutTens_d, tmpadjc1);
                  tmpfra4_2 = adimat_pop1;
                  for d=fliplr(1 : tmpfra4_2)
                     OutTens_d = adimat_pop1;
                     a_tmpca1 = adimat_adjsum(a_tmpca1, adimat_adjred(tmpca1, a_OutTens_d));
                     tmpsa1 = a_OutTens_d;
                     a_OutTens_d = a_zeros1(OutTens_d);
                     a_OutTens_d = adimat_adjsum(a_OutTens_d, adimat_adjred(OutTens_d, tmpsa1));
                     tmpca1 = adimat_pop1;
                     a_tmpca2 = adimat_adjsum(a_tmpca2, -(exp(-tmpca2) .* a_tmpca1));
                     a_tmpca1 = a_zeros1(tmpca1);
                     tmpca2 = adimat_pop1;
                     a_tmpca3 = adimat_adjsum(a_tmpca3, a_sum(a_tmpca2, tmpca3, 2));
                     a_tmpca2 = a_zeros1(tmpca2);
                     tmpca3 = adimat_pop1;
                     a_Tmmp = adimat_adjsum(a_Tmmp, a_sum(a_tmpca3, Tmmp, 1));
                     a_tmpca3 = a_zeros1(tmpca3);
                     tmpba3 = adimat_pop1;
                     if tmpba3 == 1
                        Tmmp = adimat_pop1;
                        a_tmpca1 = adimat_adjsum(a_tmpca1, adimat_adjred(tmpca1, a_Tmmp .* tmpda3));
                        a_Tmmp = a_zeros1(Tmmp);
                        tmpca1 = adimat_pop1;
                        a_W(SelfInd{n, i, v}{2}{d}, :) = adimat_adjsum(a_W(SelfInd{n, i, v}{2}{d}, :), a_repmat(a_tmpca1, W(SelfInd{n, i, v}{2}{d}, :), tmpda2));
                        a_tmpca1 = a_zeros1(tmpca1);
                        [tmpda2 tmpda3 tmpda4] = adimat_pop;
                     else
                        Tmmp = adimat_pop1;
                        a_Tmmp = a_zeros1(Tmmp);
                     end
                  end
                  OutTens_d = adimat_pop1;
                  a_OutTens_d = a_zeros1(OutTens_d);
               end
               tmpba2 = adimat_pop1;
               if tmpba2 == 1
                  InTens = adimat_pop1;
                  a_InTens = a_zeros1(InTens);
               else
                  InTens = adimat_pop1;
                  a_tmpca1 = adimat_adjsum(a_tmpca1, -(exp(-tmpca1) .* a_InTens));
                  a_InTens = a_zeros1(InTens);
                  tmpca1 = adimat_pop1;
                  a_tmpca2 = adimat_adjsum(a_tmpca2, a_sum(a_tmpca1, tmpca2));
                  a_tmpca1 = a_zeros1(tmpca1);
                  tmpca2 = adimat_pop1;
                  a_tmpca3 = adimat_adjsum(a_tmpca3, a_sum(a_tmpca2, tmpca3));
                  a_tmpca2 = a_zeros1(tmpca2);
                  tmpca3 = adimat_pop1;
                  a_W(SelfInd{n, i, v}{1}, :) = adimat_adjsum(a_W(SelfInd{n, i, v}{1}, :), adimat_adjred(W(SelfInd{n, i, v}{1}, :), a_tmpca3 .* SelfInd{n, i, v}{3}));
                  a_tmpca3 = a_zeros1(tmpca3);
               end
               v = adimat_pop1;
            end
            [msub index_sub] = adimat_pop;
         end
         [index XRF_v] = adimat_pop;
         a_XRF_v = a_zeros1(XRF_v);
      end
      sum_Tau = adimat_pop1;
      a_sum_Tau = a_zeros1(sum_Tau);
   end
   [L W] = adimat_pop;
   tmpsa1 = a_W;
   a_W = a_zeros1(W);
   a_W = adimat_adjsum(a_W, reshape(tmpsa1, size(W)));
end
