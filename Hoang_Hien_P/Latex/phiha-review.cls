\ProvidesClass{phiha-review}[2019/03/08 Review]
\NeedsTeXFormat{LaTeX2e}

\DeclareOption*{%
  \PassOptionsToClass{\CurrentOption}{article}%
}
\ProcessOptions\relax
 \LoadClass[11pt]{article}
%\LoadClass[12pt]{article}


%%%%%%%%%%
% Colors %
%%%%%%%%%%
\RequirePackage{xcolor}
\definecolor{gray}{HTML}{808080}



%%%%%%%%%
% Fonts %
%%%%%%%%%
%\RequirePackage{fontspec}
%\setmainfont[Path = fonts/,
%    BoldFont={Ubuntu-M.ttf}, 
%    ItalicFont={Ubuntu-LI.ttf},
%    BoldItalicFont={Ubuntu-BI.ttf}
%]{Ubuntu-L.ttf}


 
%%%%%%%%%%
% Footer %
%%%%%%%%%%
\RequirePackage[parfill]{parskip}
\RequirePackage{fancyhdr}
\RequirePackage{lastpage}
\RequirePackage{hyperref}

\newcommand{\footer}[3]{%
    \fancypagestyle{plain}{
        \fancyhf{} % Clear header/footer
        \renewcommand{\headrulewidth}{0pt}
        \fancyfoot[L]{\Large {\color{gray}{#1}}}
        \fancyfoot[C]{\Large {\color{gray}\href{mailto:#2}{#2}}}
        \fancyfoot[R]{\large {\color{gray}{#3}}}
    }
    \pagestyle{plain}% Set page style to plain.
}



%%%%%%%%%%
% Header %
%%%%%%%%%%
\RequirePackage{titlesec}
\newcommand{\header}[1]{%
    \makeatletter
    \def\maketitle{%
        \par{\Large{\hfil{Paper Review}\hfil}}%
        \par
        	\begin{center}
       	   \textbf{ \Large{#1} }
        	\end{center}
			%
        \par{}%%\vspace{-0.5cm}
    }
    \makeatother
    \titlelabel{\thetitle. }
    \titlespacing\section{0pt}{8pt}{0pt}
}



%%%%%%%%%%%%%%%%
% Enumerations %
%%%%%%%%%%%%%%%%
%\RequirePackage{enumerate}
%\RequirePackage[shortlabels]{enumitem}
%\setlist[enumerate]{nosep, leftmargin=*}
