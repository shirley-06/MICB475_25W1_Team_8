# MICB 475: Data Science Reserach in Microbiology
## Meeting Agenda and Minutes

### October 29, 2025
#### Agenda
- Review proposal (Aim 2, Timeline)

#### Minutes:
- Discuss with Hans about the dataset issues
- If dataset issue isn't resolved by proposal due date:
      - Make note of this in the dataset overview section
      - Disuss up to the demux data file
- For presentation (due Nov 30), the main story should be completed

#### Actions:
- Follow up with Hans about issues with the dataset
- Make prelimary figures for next meeting (Aim 1)
  
### October 22, 2025
#### Agenda
- Troubleshoot Github - unable to push/pull any code
- QIIME2 dataset update
    - after denoising/clustering and looking through the table.qzv and rep-seqs.qzv files, we lost alot of samples
    - should we pre-process data in QIIME2 or R?
- Proposal draft review
- https://docs.google.com/document/d/1PqAiz1TGuyMDlrQcFH7Er65gwhnD7v8JqtVO9qolSq4/edit?usp=sharing

#### Minutes
- Troubleshooted Github
- Qiime2 denoising/clustering issue:
    - Try another way to denoise (different parameters), save the one we already have
    - Paired ends should not decrease this much (may be overtrimmed reverse)
    - Try trimming at 132
- Pre-process in R (easier than in qiime2)
    - Eg. combining antibiotics and control
    - Rarefy by group
- Reviewed proposal draft:
    - Title should be: Exploring the Influence of Fresh and Fermented Vegetables Consumption on Gut Microbiome Diversity and Function in Women
    - Research question: Does the consumption of fresh vegetables versus fermented vegetables modulate gut microbiota dynamics and functional profiles in women over time?
    - Experimental aims:
          - Figure of how to do comparisons (Aim 1 vs Aim 2)
    - Aims specific approach:
      - Don't need to show code
      - Just step by step on what you're going to do
    - Can cite longtiduinal analysis from papers
    - Research question needs to be worked on
          - Want to compare: if FVs better than fresh vegetables at promoting a healthy gut flora
          - Word it to comapre FVs and fresh (Eg. Promotes a healthy gut flora, better functional outcomes in women)
          - If the changes stay longer than fresh veggies
    - Intro: 1-1.5 pages
          - 3 aims should be touched
          - If any similar work has been done
      
#### Action Items
- Email Evelyn about issue with QIIME2 if trimming at 132 still doesn't work
- Proposal draft:
    - Switch Aims 2 and 3 in terms of order
    - Include beta diversity in Aim 1
    - Instead of longditudinal analysis, do DeSeq instead
    - Look at relative abunadnace and cite papers
    - Cite packages in papers (from Youtube videos)
    - See if any literature: if FV have a longer impact (See less changes after WO than fresh)
    - See if any papers have done dietary interventions and see if the changes persist (major aim) and touch on the 3 aims in Intro
    - In Intro: eg. This is what they did in FV and fresh veggies, no need to describe in detail
    - Experimental aims: no need to talk about gender differences (not enough samples to look at both males and females, can say looking at it in context of women and its important because its sex specific) --> Looking at gut health in the context of women only, not a lot of research on women only.
    - Have stats and steps of what you're going to do
    - Timeline: give yourself more time (buffer time), overlap multiple aims in mutliple steps
    - Dataset overview: Make sure each question answered


### October 15, 2025
#### Agenda
- QIIME2 dataset processing updates
- Questions about workflow for paired end sequences
- Trimming based on demux.qzv plot

#### Minutes
- Follow paired-end reads protocol on Qiime site
- Upload Qiime code onto Github figure out how to push and pull
- Analyzing before and after eating fresh and fermented veggies
- Washout period is important to see if the effects on micromebiodiversity lasts
- Compare points 1,3,5 to see before one and after
- Research question: Does fresh vs fermentated vegetation alter microbiome composition and function in women?
- Can do core analysis on microbiome and do further research in literature to see if they are beneficial
- Combine control and antibiotics into one group, vs constipation, make sure to explain logic
- Alpha and beta diversity, including washouts, using pie crust
- Whoever does more of the coding should do less of the proposal writing and R
- When making presentation, make the timepoints clear, if they end up being irrelevant, label them n.s.
- Proposal: concise introduction, ONE hypothesis: timepoints affect microbiome diversity differently across different conditions?

#### Action Items
- - For Proposal:
- Draw analysis plans for longitudinal research and explain why
- Background with citations, hypothesis based on existing literature
- Dataset overview: Finish data processing up to filtering
- Figures with subbings? using a pie crust
- Have at least bullet points of proposal by next meeting
- Cite the R package sources


### October 8, 2025
From chad: 
look at the analyses done in this longitudinal study
https://pmc.ncbi.nlm.nih.gov/articles/PMC10832965/#s0004

#### Agenda
- Confirm project 2 dataset
- Discuss project 2 timeline

#### Minutes

Project 2 dataset: [Impact of fresh and fermented vegetable consumption on gut microbiota and body composition: insights from diverse data analysis approaches](https://www.frontiersin.org/journals/nutrition/articles/10.3389/fnut.2025.1623710/full#supplementary-material)


- Reviewed study design and microbiome analysis methodology from the paper 
    - Paper mainly looked at relative abundance diversity
    - Lots of individual data present
    - Identified gaps for deeper exploration
- Brainstormed areas of focus for project 2
    - Perform a more in-depth analysis than the original study
    - Potential analytical directions:  
        - Core analysis
        - Indicator analysis
        - Functional analysis (performed using picrust2 - see module 19)
        - Paired analysis
        - Longitudinal analysis
    - Use strident cutoffs, controlling for sex, antibiotics, constipation
    - Compare pre/post dietary intervention
    - Investigate changes during washout periods - do the communities persist?
- Project 2 proposal is due Oct 26
      - QIIME2 workflow needs to be completed for the proposal

#### Action Items
- Look through the dataset and begin QIIME2 workflow (Dr. Sun will notify us it has been uploaded to the server)
- Begin background research for the proposal
- Review all the suggested analysis methods
    -  Review longitudinal analysis from the paper Chad shared 




### October 1, 2025
#### Minutes
- For each meeting:
    - Prepare an agenda
        - ensure to complete the night before
    - Assign someone to take meeting minutes
 - Discussed interest in medical microbiome datasets
     - Briefly looked at a few datasets from NASA data repository
  
#### Action Items
- Explore datasets and brainstorm potential research question
    - Compile new datasets and send them to Dr. Sun Friday EOD
- Add Chad to the repository
