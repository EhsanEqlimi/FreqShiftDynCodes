                                                                    %% =========================================
% MATLAB → GitHub CLEAN RESTART (AFTER .git DELETION)
% FIXES: broken repo + large file prevention + push conflict
%% =========================================

%% 1. Go to project folder
cd('/Users/ehsaneqlimi/FreqShiftDynCodes')
% ⬆️ change if needed

%% 2. Check files
!ls -a

%% =========================================
% 3. INIT NEW GIT REPO
%% =========================================
!git init

%% 4. Set branch name
!git branch -M main

%% =========================================
% 5. IMPORTANT: prevent large file problem
%% =========================================
!echo "*.mat" >> .gitignore

%% =========================================
% 6. Add files
%% =========================================
!git add .

%% 7. First commit (safe check)
%% =========================================
!git commit -m "Clean restart repo after .git deletion" || echo "Nothing to commit"

%% =========================================
% 8. Add GitHub remote (SSH)
%% =========================================
!git remote remove origin
!git remote add origin git@github.com:EhsanEqlimi/FreqShiftDynCodes.git

%% =========================================
% 9. PUSH (FIXED)
% This solves: "remote contains work you do not have locally"
%% =========================================
!git push -u origin main --force