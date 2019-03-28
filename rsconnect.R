# rsconnect::setAccountInfo(
#   name='catholic', 
#   token='35A4CF5EFA0831AA195176F3F0C2A40D', 
#   secret='Onpx1AJJQ5HXlxqdTzQ/O1jODgLF0075KzGV9flR'
# )

# https://github.com/pipetcpt/SWD-meeting/wiki/Deploying-a-Shiny-app

rsconnect::setAccountInfo(
  name='pipet',
  token='어쩌고 저쩌고',
  secret='어쩌고 저쩌고'
)
rsconnect::deployApp(account = "pipet", appName = "mydrug")