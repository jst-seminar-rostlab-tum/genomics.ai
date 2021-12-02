import React from 'react';
import styles from './about.module.css';
import { Typography ,Card,CardMedia,CardContent ,Grid, Button,Link,IconButton, Divider} from '@mui/material';
import NavBar from '../../../NavBar/NavBar';
import Footer from '../../Footer/Footer';
import FacebookIcon from '@mui/icons-material/Facebook';
import TwitterIcon from '@mui/icons-material/Twitter';
import GitHubIcon from '@mui/icons-material/GitHub';
import LinkedInIcon from '@mui/icons-material/LinkedIn';



function Namecard(props) {
  return(
    <Card sx={{ maxWidth: 500,
      display: 'flex',
      paddingLeft:'20px',
      paddingright:'20px',
      boxShadow: "none"
      }}>
    <CardMedia
      component="img"
      image={props.img}
      alt={props.name}
      sx={{
        borderRadius: '50%',
        height: '120px',
        width: '120px',
        margin: '28px',
        objectFit: 'cover',
        border:'4px solid #ff',
        display: 'flex',
      }}
    />

    <CardContent style={{ maxWidth: "200px"}}>
      <Typography variant = "h5">{props.name}</Typography>
      <Typography variant="body2" color="text.secondary">{props.role}</Typography>
      <br/>
      <Typography variant="body2" 
      style={{ maxWidth: '200px',
          alignItems: 'center',
          wordWrap: 'break-word'
          }}
      >{props.dscp}</Typography>

      <Grid container spacing={0}
        sx = {{ 
          marginTop: '0px',
          justifyContent: 'center',
        }}
      >
        
      {
        props.socialFB != ''
        ? FbItem(props.socialFB)
        : <Grid item/>
      }

      {
        props.socialGithub != ''
        ? GithubItem(props.socialGithub)
        : <Grid item/>
      }
      
      {
        props.socialLinkedIn != ''
        ? LinkedInItem(props.socialLinkedIn)
        : <Grid item/>
      }

      {
        props.socialTwitter != ''
        ? TwitterItem(props.socialTwitter)
        : <Grid item/>
      }

      </Grid>

    </CardContent>
  </Card>

  );
}

function FbItem(fblink){
    return (
      <Grid item>
        <IconButton href={fblink}>
          <FacebookIcon style={{fill : '#3B8EED'}}/>
        </IconButton>
        </Grid>
    )
}

function GithubItem(ghlink){
  return (
    <Grid item>
      <IconButton href={ghlink}>
        <GitHubIcon style={{fill : '#171B22'}}/>
      </IconButton>
      </Grid>
  )
}

function LinkedInItem(lilink){
  return (
    <Grid item>
      <IconButton href={lilink}>
        <LinkedInIcon style={{fill : '#2D69BF'}}/>
      </IconButton>
      </Grid>
  )
}

function TwitterItem(ttlink){
  return (
    <Grid item>
      <IconButton href={ttlink}>
        <TwitterIcon style={{fill : '#66AFEC'}}/>
      </IconButton>
      </Grid>
  )
}

const SubteamSection = ((props) => {
  const cardWidth = 3.5;
  const teamlead = props.data.shift();
  console.log(teamlead.img);
  return(
    <div>

    <Grid
      container
      spacing={1}
      direction="row"
      justifyContent="center"
      // justifyContent="space-between"
      alignItems="center"
      justify="center"
      style={{ minHeight: '20vh', maxWidth: '500vh' }}
      sx={{paddingTop:'50px'}}
      >
        <Grid item xs={cardWidth}>
          <Namecard 
          name = {teamlead.name}
          role={teamlead.role}
          img={teamlead.img}
          dscp={teamlead.dscp}
          socialFB={teamlead.socialFB}
          socialLinkedIn={teamlead.socialLinkedIn}
          socialGithub={teamlead.socialGithub}
          socialTwitter={teamlead.socialTwitter}
          />
        </Grid>  

    </Grid>
    
    <Grid
    container
    spacing={2}
    direction="row"
    justifyContent="center"
    // justifyContent="space-between"
    alignItems="center"
    justify="center"
    style={{ minHeight: '20vh', maxWidth: '500vh' }}
    sx={{paddingTop:'50px'}}
    >
    
      {
        props.data.map(elem =>(
          <Grid item xs={cardWidth}>
            <Namecard
            name = {elem.name}
            role={elem.role}
            img={elem.img}
            dscp={elem.dscp}
            socialFB={elem.socialFB}
            socialLinkedIn={elem.socialLinkedIn}
            socialGithub={elem.socialGithub}
            socialTwitter={elem.socialTwitter}
            />

          </Grid>
        ))
        
        }

    </Grid>
    </div>
  )

}
)

function About() {
  const cardWidth = 3.5;
  const frontEndData = [
    {name: 'Amin Ben Saad', role : 'Frontend Lead & Developer' , img : 'memberphotos/Frontend_Amin Ben Saad.png' , dscp: 'B.Sc. Informatics at TUM ' ,socialFB:'',socialGithub:'',socialTwitter:'',socialLinkedIn:'https://www.linkedin.com/in/amin-ben-saad369/?originalSubdomain=de'},
    {name: 'Dominik Jámbor', role : 'Frontend Developer' , img : 'memberphotos/Dominik.jpg' , dscp: 'Full-stack web developer and M.Sc. student at ELTE Budapest. ' ,socialFB:'',socialGithub:'https://github.com/dominikjambor',socialTwitter:'',socialLinkedIn:''},
    {name: 'Ronald Skorobogat', role : 'Frontend Developer' , img : 'memberphotos/Frontend-RonaldSkorobogat.jpg' , dscp: 'B.Sc. Informatics at TUM ' ,socialFB:'',socialGithub:'https://github.com/ronskoro',socialTwitter:'',socialLinkedIn:'https://www.linkedin.com/in/ronald-skorobogat-a1348b1a3/'},
    {name: 'Liudongnan Yang', role : 'Frontend Developer' , img : 'memberphotos/Yang.jpg' , dscp: " Master's student in Informatics at TUM" ,socialFB:'',socialGithub:'https://github.com/yldn',socialTwitter:'',socialLinkedIn:'https://www.linkedin.com/in/liudongnan-yang-683a3b163/'}
  ];

  const backEndData = [
    {name: 'Lars Frölich', role : 'Backend Lead & Developer' , img : '' , dscp: 'B.Sc. Informatics at TUM ' ,socialFB:'',socialGithub:'',socialTwitter:'',socialLinkedIn:''},
    {name: 'Matthias Michailow' , role : 'Backend Developer' , img : 'memberphotos/Backend_Matthias Michailow.jpg' , dscp: "Master's student in Informatics at TUM " ,socialFB:'',socialGithub:'',socialTwitter:'',socialLinkedIn:''},
    {name: 'Yi Rui Cui', role : 'Backend Developer' , img : '' , dscp: 'B.Sc. Informatics: Games Engineering at TUM' ,socialFB:'',socialGithub:'',socialTwitter:'',socialLinkedIn:''},
    {name: 'Manuel Brandstetter', role : 'Backend Developer' , img : '' , dscp: " Full-stack web developer and Master's student in Informatics at TUM" ,socialFB:'',socialGithub:'',socialTwitter:'',socialLinkedIn:''}
  ];

  const visData = [
    {name: 'Vivian Haller', role : 'Visualisation Team lead, d3.js developer' , img : 'https://kiranworkspace.com/demo/projects/code-snippets/team/our-team2/images/users/user1.jpg' , dscp: 'M.Sc. Math/CS at TUM' ,socialFB:'',socialGithub:'',socialTwitter:'',socialLinkedIn:''},
    {name: 'Yashica Saun' , role : 'Visualisation/d3.js Developer' , img : 'https://kiranworkspace.com/demo/projects/code-snippets/team/our-team2/images/users/user1.jpg' , dscp: "Master's student in Informatics at TUM" ,socialFB:'',socialGithub:'',socialTwitter:'',socialLinkedIn:''},
    {name: 'Carl Neumann', role : 'Visualisation/d3.js Developer' , img : 'https://kiranworkspace.com/demo/projects/code-snippets/team/our-team2/images/users/user1.jpg' , dscp: 'B.Sc. Information Systems at TUM ' ,socialFB:'',socialGithub:'',socialTwitter:'',socialLinkedIn:''},
    {name: 'Mohan Prabhakaran', role : 'Visualisation/d3.js Developer' , img : 'https://kiranworkspace.com/demo/projects/code-snippets/team/our-team2/images/users/user1.jpg' , dscp: "Master's student in Informatics at TUM" ,socialFB:'',socialGithub:'',socialTwitter:'',socialLinkedIn:''}
  ];

  return (
    <div className={styles.headerContainer}>
      <NavBar />

      <Typography sx={{ fontWeight: '400', fontSize: '24px' }}>Team</Typography>
      
      
{/* For Guy */}
          <Divider variant="middle" textAlign = "center" sx={{padding: '20px'}}></Divider>

          <Grid
            container
            spacing={0}
            direction="column"
            alignItems="center"
            justifyContent="center"
            // style={{ minHeight: '100vh' }}
          >
            <Grid item xs={3}>
              <Namecard 
              name="Guy Yachdav" 
              role="Supervisor & Initiator" 
              img='https://scholar.googleusercontent.com/citations?view_op=medium_photo&user=UoUkGhUAAAAJ&citpid=2' 
              dscp='Good looking guy with a lot of experience'
              socialFB=''
              socialLinkedIn='https://www.linkedin.com/in/gyachdav/?originalSubdomain=il'
              socialGithub=''
              socialTwitter=''
              />
            </Grid>

          </Grid>

          <text>GeneCruncher was developed by the following team</text>
          <br />
          <Divider variant="middle" textAlign = "left" sx={{padding: '20px'}}>Frontend Team</Divider>
          
          <SubteamSection data={frontEndData} />


          <Divider variant="middle" textAlign = "left" sx={{padding: '20px'}}>Backend Team</Divider>

          <SubteamSection data={backEndData} />


          <Divider variant="middle" textAlign = "left" sx={{padding: '20px'}}>Visualization Team</Divider>

          <SubteamSection data={visData} />

      <Footer />
    </div>
  );
}

export default About;