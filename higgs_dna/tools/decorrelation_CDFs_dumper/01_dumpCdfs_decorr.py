{"payload":{"allShortcutsEnabled":false,"fileTree":{"":{"items":[{"name":"01_dumpCdfs_decorr.py","path":"01_dumpCdfs_decorr.py","contentType":"file"},{"name":"02_decorrelate.py","path":"02_decorrelate.py","contentType":"file"},{"name":"README.md","path":"README.md","contentType":"file"},{"name":"plot_utils_smear_term.py","path":"plot_utils_smear_term.py","contentType":"file"},{"name":"ploter.py","path":"ploter.py","contentType":"file"}],"totalCount":5}},"fileTreeProcessingTime":1.9829829999999997,"foldersToFetch":[],"repo":{"id":734388427,"defaultBranch":"main","name":"CDF_decorrelation","ownerLogin":"CaioDaumann","currentUserCanPush":false,"isFork":false,"isEmpty":false,"createdAt":"2023-12-21T14:57:45.000Z","ownerAvatar":"https://avatars.githubusercontent.com/u/146199123?v=4","public":true,"private":false,"isOrgOwned":false},"symbolsExpanded":false,"treeExpanded":true,"refInfo":{"name":"main","listCacheKey":"v0:1703170774.0","canEdit":false,"refType":"branch","currentOid":"6f8fc1a501c999d0f4e900c60ae9a7b4522aaafd"},"path":"01_dumpCdfs_decorr.py","currentUser":null,"blob":{"rawLines":["# import root_pandas","import numpy as np","import pandas as pd","from decorrelator import cdfCalc","import argparse","from coffea.nanoevents import NanoEventsFactory, BaseSchema","import os","","def main(options):","","    df = pd.DataFrame()","","    import glob","    #files = glob.glob( str(options.infile) + \"*.parquet\")","    #data = [pd.read_parquet(f) for f in files]","    #data = pd.read_parquet(\"/net/scratch_cms3a/daumann/massresdecorrhiggsdna/big_bkg/Diphoton.parquet\")","    #events= pd.concat(data,ignore_index=True)","","    files  = glob.glob( \"/net/scratch_cms3a/daumann/HiggsDNA/diphoton_samples_for_CDFs/Diphoton_2022_postEE/nominal/*.parquet\" )","    data   = [pd.read_parquet(f) for f in files]","    events = pd.concat(data,ignore_index=True)","","    df[\"sigma_m_over_m\"] = events.sigma_m_over_m_Smeared.to_numpy()","    df[\"mass\"]   = events.mass.to_numpy()","    df[\"weight\"] = events.weight.to_numpy()    ","","    print(f\"INFO: found {len(events)} events\")","","    # clarify with anyone what this \"sigmarv\" business here is about!","    df[\"sigma_m_over_m\"] = events.sigma_m_over_m_Smeared.to_numpy()","","    df[\"mass\"]   = events.mass.to_numpy()","    df[\"weight\"] = events.weight.to_numpy()","","    # df = root_pandas.read_root(options.infile, options.tree, columns=['sigmarv', 'recoMass', 'weight'])","    calc = cdfCalc(df, 'sigma_m_over_m', 'mass', np.linspace(100, 180, 161)) # 161 bins initially","    calc.dumpCdfs(options.cdfsFile)","","","if __name__ == \"__main__\":","","    parser = argparse.ArgumentParser()","    requiredArgs = parser.add_argument_group()","    #requiredArgs.add_argument('-i', '--infile', action='store', type=str, required=True)","    requiredArgs.add_argument('-i', '--infile', action='store', type=str, required=True)","    requiredArgs.add_argument('-c', '--cdfsFile', action='store', type=str, required=True)","    requiredArgs.add_argument('-t', '--tree', action='store', type=str, required=True)","    options = parser.parse_args()","    main(options)"],"stylingDirectives":[[{"start":0,"end":20,"cssClass":"pl-c"}],[{"start":0,"end":6,"cssClass":"pl-k"},{"start":7,"end":12,"cssClass":"pl-s1"},{"start":13,"end":15,"cssClass":"pl-k"},{"start":16,"end":18,"cssClass":"pl-s1"}],[{"start":0,"end":6,"cssClass":"pl-k"},{"start":7,"end":13,"cssClass":"pl-s1"},{"start":14,"end":16,"cssClass":"pl-k"},{"start":17,"end":19,"cssClass":"pl-s1"}],[{"start":0,"end":4,"cssClass":"pl-k"},{"start":5,"end":17,"cssClass":"pl-s1"},{"start":18,"end":24,"cssClass":"pl-k"},{"start":25,"end":32,"cssClass":"pl-s1"}],[{"start":0,"end":6,"cssClass":"pl-k"},{"start":7,"end":15,"cssClass":"pl-s1"}],[{"start":0,"end":4,"cssClass":"pl-k"},{"start":5,"end":11,"cssClass":"pl-s1"},{"start":12,"end":22,"cssClass":"pl-s1"},{"start":23,"end":29,"cssClass":"pl-k"},{"start":30,"end":47,"cssClass":"pl-v"},{"start":49,"end":59,"cssClass":"pl-v"}],[{"start":0,"end":6,"cssClass":"pl-k"},{"start":7,"end":9,"cssClass":"pl-s1"}],[],[{"start":0,"end":3,"cssClass":"pl-k"},{"start":4,"end":8,"cssClass":"pl-en"},{"start":9,"end":16,"cssClass":"pl-s1"}],[],[{"start":4,"end":6,"cssClass":"pl-s1"},{"start":7,"end":8,"cssClass":"pl-c1"},{"start":9,"end":11,"cssClass":"pl-s1"},{"start":12,"end":21,"cssClass":"pl-v"}],[],[{"start":4,"end":10,"cssClass":"pl-k"},{"start":11,"end":15,"cssClass":"pl-s1"}],[{"start":4,"end":58,"cssClass":"pl-c"}],[{"start":4,"end":47,"cssClass":"pl-c"}],[{"start":4,"end":104,"cssClass":"pl-c"}],[{"start":4,"end":46,"cssClass":"pl-c"}],[],[{"start":4,"end":9,"cssClass":"pl-s1"},{"start":11,"end":12,"cssClass":"pl-c1"},{"start":13,"end":17,"cssClass":"pl-s1"},{"start":18,"end":22,"cssClass":"pl-en"},{"start":24,"end":126,"cssClass":"pl-s"}],[{"start":4,"end":8,"cssClass":"pl-s1"},{"start":11,"end":12,"cssClass":"pl-c1"},{"start":14,"end":16,"cssClass":"pl-s1"},{"start":17,"end":29,"cssClass":"pl-en"},{"start":30,"end":31,"cssClass":"pl-s1"},{"start":33,"end":36,"cssClass":"pl-k"},{"start":37,"end":38,"cssClass":"pl-s1"},{"start":39,"end":41,"cssClass":"pl-c1"},{"start":42,"end":47,"cssClass":"pl-s1"}],[{"start":4,"end":10,"cssClass":"pl-s1"},{"start":11,"end":12,"cssClass":"pl-c1"},{"start":13,"end":15,"cssClass":"pl-s1"},{"start":16,"end":22,"cssClass":"pl-en"},{"start":23,"end":27,"cssClass":"pl-s1"},{"start":28,"end":40,"cssClass":"pl-s1"},{"start":40,"end":41,"cssClass":"pl-c1"},{"start":41,"end":45,"cssClass":"pl-c1"}],[],[{"start":4,"end":6,"cssClass":"pl-s1"},{"start":7,"end":23,"cssClass":"pl-s"},{"start":25,"end":26,"cssClass":"pl-c1"},{"start":27,"end":33,"cssClass":"pl-s1"},{"start":34,"end":56,"cssClass":"pl-s1"},{"start":57,"end":65,"cssClass":"pl-en"}],[{"start":4,"end":6,"cssClass":"pl-s1"},{"start":7,"end":13,"cssClass":"pl-s"},{"start":17,"end":18,"cssClass":"pl-c1"},{"start":19,"end":25,"cssClass":"pl-s1"},{"start":26,"end":30,"cssClass":"pl-s1"},{"start":31,"end":39,"cssClass":"pl-en"}],[{"start":4,"end":6,"cssClass":"pl-s1"},{"start":7,"end":15,"cssClass":"pl-s"},{"start":17,"end":18,"cssClass":"pl-c1"},{"start":19,"end":25,"cssClass":"pl-s1"},{"start":26,"end":32,"cssClass":"pl-s1"},{"start":33,"end":41,"cssClass":"pl-en"}],[],[{"start":4,"end":9,"cssClass":"pl-en"},{"start":10,"end":45,"cssClass":"pl-s"},{"start":24,"end":37,"cssClass":"pl-s1"},{"start":24,"end":25,"cssClass":"pl-kos"},{"start":25,"end":28,"cssClass":"pl-en"},{"start":29,"end":35,"cssClass":"pl-s1"},{"start":36,"end":37,"cssClass":"pl-kos"}],[],[{"start":4,"end":69,"cssClass":"pl-c"}],[{"start":4,"end":6,"cssClass":"pl-s1"},{"start":7,"end":23,"cssClass":"pl-s"},{"start":25,"end":26,"cssClass":"pl-c1"},{"start":27,"end":33,"cssClass":"pl-s1"},{"start":34,"end":56,"cssClass":"pl-s1"},{"start":57,"end":65,"cssClass":"pl-en"}],[],[{"start":4,"end":6,"cssClass":"pl-s1"},{"start":7,"end":13,"cssClass":"pl-s"},{"start":17,"end":18,"cssClass":"pl-c1"},{"start":19,"end":25,"cssClass":"pl-s1"},{"start":26,"end":30,"cssClass":"pl-s1"},{"start":31,"end":39,"cssClass":"pl-en"}],[{"start":4,"end":6,"cssClass":"pl-s1"},{"start":7,"end":15,"cssClass":"pl-s"},{"start":17,"end":18,"cssClass":"pl-c1"},{"start":19,"end":25,"cssClass":"pl-s1"},{"start":26,"end":32,"cssClass":"pl-s1"},{"start":33,"end":41,"cssClass":"pl-en"}],[],[{"start":4,"end":105,"cssClass":"pl-c"}],[{"start":4,"end":8,"cssClass":"pl-s1"},{"start":9,"end":10,"cssClass":"pl-c1"},{"start":11,"end":18,"cssClass":"pl-en"},{"start":19,"end":21,"cssClass":"pl-s1"},{"start":23,"end":39,"cssClass":"pl-s"},{"start":41,"end":47,"cssClass":"pl-s"},{"start":49,"end":51,"cssClass":"pl-s1"},{"start":52,"end":60,"cssClass":"pl-en"},{"start":61,"end":64,"cssClass":"pl-c1"},{"start":66,"end":69,"cssClass":"pl-c1"},{"start":71,"end":74,"cssClass":"pl-c1"},{"start":77,"end":97,"cssClass":"pl-c"}],[{"start":4,"end":8,"cssClass":"pl-s1"},{"start":9,"end":17,"cssClass":"pl-en"},{"start":18,"end":25,"cssClass":"pl-s1"},{"start":26,"end":34,"cssClass":"pl-s1"}],[],[],[{"start":0,"end":2,"cssClass":"pl-k"},{"start":3,"end":11,"cssClass":"pl-s1"},{"start":12,"end":14,"cssClass":"pl-c1"},{"start":15,"end":25,"cssClass":"pl-s"}],[],[{"start":4,"end":10,"cssClass":"pl-s1"},{"start":11,"end":12,"cssClass":"pl-c1"},{"start":13,"end":21,"cssClass":"pl-s1"},{"start":22,"end":36,"cssClass":"pl-v"}],[{"start":4,"end":16,"cssClass":"pl-s1"},{"start":17,"end":18,"cssClass":"pl-c1"},{"start":19,"end":25,"cssClass":"pl-s1"},{"start":26,"end":44,"cssClass":"pl-en"}],[{"start":4,"end":89,"cssClass":"pl-c"}],[{"start":4,"end":16,"cssClass":"pl-s1"},{"start":17,"end":29,"cssClass":"pl-en"},{"start":30,"end":34,"cssClass":"pl-s"},{"start":36,"end":46,"cssClass":"pl-s"},{"start":48,"end":54,"cssClass":"pl-s1"},{"start":54,"end":55,"cssClass":"pl-c1"},{"start":55,"end":62,"cssClass":"pl-s"},{"start":64,"end":68,"cssClass":"pl-s1"},{"start":68,"end":69,"cssClass":"pl-c1"},{"start":69,"end":72,"cssClass":"pl-s1"},{"start":74,"end":82,"cssClass":"pl-s1"},{"start":82,"end":83,"cssClass":"pl-c1"},{"start":83,"end":87,"cssClass":"pl-c1"}],[{"start":4,"end":16,"cssClass":"pl-s1"},{"start":17,"end":29,"cssClass":"pl-en"},{"start":30,"end":34,"cssClass":"pl-s"},{"start":36,"end":48,"cssClass":"pl-s"},{"start":50,"end":56,"cssClass":"pl-s1"},{"start":56,"end":57,"cssClass":"pl-c1"},{"start":57,"end":64,"cssClass":"pl-s"},{"start":66,"end":70,"cssClass":"pl-s1"},{"start":70,"end":71,"cssClass":"pl-c1"},{"start":71,"end":74,"cssClass":"pl-s1"},{"start":76,"end":84,"cssClass":"pl-s1"},{"start":84,"end":85,"cssClass":"pl-c1"},{"start":85,"end":89,"cssClass":"pl-c1"}],[{"start":4,"end":16,"cssClass":"pl-s1"},{"start":17,"end":29,"cssClass":"pl-en"},{"start":30,"end":34,"cssClass":"pl-s"},{"start":36,"end":44,"cssClass":"pl-s"},{"start":46,"end":52,"cssClass":"pl-s1"},{"start":52,"end":53,"cssClass":"pl-c1"},{"start":53,"end":60,"cssClass":"pl-s"},{"start":62,"end":66,"cssClass":"pl-s1"},{"start":66,"end":67,"cssClass":"pl-c1"},{"start":67,"end":70,"cssClass":"pl-s1"},{"start":72,"end":80,"cssClass":"pl-s1"},{"start":80,"end":81,"cssClass":"pl-c1"},{"start":81,"end":85,"cssClass":"pl-c1"}],[{"start":4,"end":11,"cssClass":"pl-s1"},{"start":12,"end":13,"cssClass":"pl-c1"},{"start":14,"end":20,"cssClass":"pl-s1"},{"start":21,"end":31,"cssClass":"pl-en"}],[{"start":4,"end":8,"cssClass":"pl-en"},{"start":9,"end":16,"cssClass":"pl-s1"}]],"csv":null,"csvError":null,"dependabotInfo":{"showConfigurationBanner":false,"configFilePath":null,"networkDependabotPath":"/CaioDaumann/CDF_decorrelation/network/updates","dismissConfigurationNoticePath":"/settings/dismiss-notice/dependabot_configuration_notice","configurationNoticeDismissed":null,"repoAlertsPath":"/CaioDaumann/CDF_decorrelation/security/dependabot","repoSecurityAndAnalysisPath":"/CaioDaumann/CDF_decorrelation/settings/security_analysis","repoOwnerIsOrg":false,"currentUserCanAdminRepo":false},"displayName":"01_dumpCdfs_decorr.py","displayUrl":"https://github.com/CaioDaumann/CDF_decorrelation/blob/main/01_dumpCdfs_decorr.py?raw=true","headerInfo":{"blobSize":"1.88 KB","deleteInfo":{"deleteTooltip":"You must be signed in to make or propose changes"},"editInfo":{"editTooltip":"You must be signed in to make or propose changes"},"ghDesktopPath":"https://desktop.github.com","gitLfsPath":null,"onBranch":true,"shortPath":"0599b65","siteNavLoginPath":"/login?return_to=https%3A%2F%2Fgithub.com%2FCaioDaumann%2FCDF_decorrelation%2Fblob%2Fmain%2F01_dumpCdfs_decorr.py","isCSV":false,"isRichtext":false,"toc":null,"lineInfo":{"truncatedLoc":"49","truncatedSloc":"37"},"mode":"file"},"image":false,"isCodeownersFile":null,"isPlain":false,"isValidLegacyIssueTemplate":false,"issueTemplateHelpUrl":"https://docs.github.com/articles/about-issue-and-pull-request-templates","issueTemplate":null,"discussionTemplate":null,"language":"Python","languageID":303,"large":false,"loggedIn":false,"newDiscussionPath":"/CaioDaumann/CDF_decorrelation/discussions/new","newIssuePath":"/CaioDaumann/CDF_decorrelation/issues/new","planSupportInfo":{"repoIsFork":null,"repoOwnedByCurrentUser":null,"requestFullPath":"/CaioDaumann/CDF_decorrelation/blob/main/01_dumpCdfs_decorr.py","showFreeOrgGatedFeatureMessage":null,"showPlanSupportBanner":null,"upgradeDataAttributes":null,"upgradePath":null},"publishBannersInfo":{"dismissActionNoticePath":"/settings/dismiss-notice/publish_action_from_dockerfile","releasePath":"/CaioDaumann/CDF_decorrelation/releases/new?marketplace=true","showPublishActionBanner":false},"rawBlobUrl":"https://github.com/CaioDaumann/CDF_decorrelation/raw/main/01_dumpCdfs_decorr.py","renderImageOrRaw":false,"richText":null,"renderedFileInfo":null,"shortPath":null,"symbolsEnabled":true,"tabSize":8,"topBannersInfo":{"overridingGlobalFundingFile":false,"globalPreferredFundingPath":null,"repoOwner":"CaioDaumann","repoName":"CDF_decorrelation","showInvalidCitationWarning":false,"citationHelpUrl":"https://docs.github.com/github/creating-cloning-and-archiving-repositories/creating-a-repository-on-github/about-citation-files","showDependabotConfigurationBanner":false,"actionsOnboardingTip":null},"truncated":false,"viewable":true,"workflowRedirectUrl":null,"symbols":{"timed_out":false,"not_analyzed":false,"symbols":[{"name":"main","kind":"function","ident_start":184,"ident_end":188,"extent_start":180,"extent_end":1399,"fully_qualified_name":"main","ident_utf16":{"start":{"line_number":8,"utf16_col":4},"end":{"line_number":8,"utf16_col":8}},"extent_utf16":{"start":{"line_number":8,"utf16_col":0},"end":{"line_number":36,"utf16_col":35}}}]}},"copilotInfo":null,"copilotAccessAllowed":false,"csrf_tokens":{"/CaioDaumann/CDF_decorrelation/branches":{"post":"qM3f0Cvb2Lx5M6viDiKrI7euG_LmPQJoNyAtRiYm4zzIljJIP2nfctQ6rM-9JyGtXzMVljKsRDz_N1zcsgar-A"},"/repos/preferences":{"post":"M6mTu1sg2Bin9bJ97JkvMGO0LAbym0_Iuf1m8pA9b4mrWh5bOg9-7wQ_mkAcXor3fNG8FA7TChbNbud8UWjb0g"}}},"title":"CDF_decorrelation/01_dumpCdfs_decorr.py at main · CaioDaumann/CDF_decorrelation"}