#clear github actions caches
#list cache IDs
ghjson = system(paste0('curl -H "Accept: application/vnd.github+json" -H "Authorization: Bearer ',Sys.getenv('GITHUB_PAT'),'" https://api.github.com/repos/pboesu/moultmcmc/actions/caches'))

#delete cache IDs (adjust loop index from above)
for (id in 319:324) {
ghjson = system(paste0('curl -X DELETE -H "Accept: application/vnd.github+json" -H "Authorization: Bearer ',Sys.getenv('GITHUB_PAT'),'" https://api.github.com/repos/pboesu/moultmcmc/actions/caches/',id))
}
