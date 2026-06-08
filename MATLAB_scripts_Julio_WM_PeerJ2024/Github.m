%% ==============================
%  MATLAB → GitHub (FIRST PUSH SETUP)
%  Requires: GitHub repo already created
%% ==============================

%% 1. Go to your project folder
cd('/Users/ehsaneqlimi/FreqShiftDynCodes')
% ⬆️ change if needed

%% 2. Check git
!git --version

%% 3. Initialize repository (safe even if already exists)
!git init

%% 4. Add all files
!git add .

%% 5. Commit
!git commit -m "Initial commit"

%% 6. Remove any old remote (safe)
!git remote remove origin

%% 7. Add NEW GitHub SSH remote
!git remote add origin git@github.com:EhsanEqlimi/FreqShiftDynCodes.git

%% 8. Ensure branch is main
!git branch -M main

%% 9. Test SSH (optional but useful)
!ssh -T git@github.com

%% 10. Push to GitHub
!git push -u origin main