import React from 'react';
import styles from './about.module.css';
import { Typography ,Card,CardMedia,CardContent ,Grid, Button,Link,IconButton} from '@mui/material';
import NavBar from '../../../NavBar/NavBar';
import FacebookIcon from '@mui/icons-material/Facebook';
import TwitterIcon from '@mui/icons-material/Twitter';
import GitHubIcon from '@mui/icons-material/GitHub';
import LinkedInIcon from '@mui/icons-material/LinkedIn';

function MediaCards(props) {
    const { classes } = props;
    const cardWidth = 3.5;
    return(
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
    
        <Grid item xs={cardWidth}>
          <Namecard 
          name="Paella dish" 
          role="CEO, Co-founder" 
          img='https://kiranworkspace.com/demo/projects/code-snippets/team/our-team2/images/users/user1.jpg' 
          dscp='BlablaggQWQEQWadweqweqwdqwd123 ajhwlkellkcnl adq'
          socialFB='#!'
          socialLinkedIn='#!'
          socialGithub='#!'
          socialTwitter='#!'
          />
        </Grid>  
        <Grid item xs={cardWidth}>
          <Namecard 
          name="Paella dish" 
          role="CEO, Co-founder" 
          img='https://kiranworkspace.com/demo/projects/code-snippets/team/our-team2/images/users/user1.jpg' 
          dscp='BlablaggQWQEQWadweqweqwdqwd123 ajhwlkellkcnl adq'
          socialFB='#!'
          socialLinkedIn='#!'
          socialGithub='#!'
          socialTwitter='#!'
          />
        </Grid>  
        <Grid item xs={cardWidth}>
          <Namecard 
          name="Paella dish" 
          role="CEO, Co-founder" 
          img='https://kiranworkspace.com/demo/projects/code-snippets/team/our-team2/images/users/user1.jpg' 
          dscp='BlablaggQWQEQWadweqweqwdqwd123 ajhwlkellkcnl adq'
          socialFB='#!'
          socialLinkedIn='#!'
          socialGithub='#!'
          socialTwitter='#!'
          />
        </Grid>  
        <Grid item xs={cardWidth}>
          <Namecard 
          name="Paella dish" 
          role="CEO, Co-founder" 
          img='https://kiranworkspace.com/demo/projects/code-snippets/team/our-team2/images/users/user1.jpg' 
          dscp='BlablaggQWQEQWadweqweqwdqwd123 ajhwlkellkcnl adq'
          socialFB='#!'
          socialLinkedIn='#!'
          socialGithub='#!'
          socialTwitter='#!'
          />
        </Grid>  
        <Grid item xs={cardWidth}>
          <Namecard 
          name="Paella dish" 
          role="CEO, Co-founder" 
          img='https://kiranworkspace.com/demo/projects/code-snippets/team/our-team2/images/users/user1.jpg' 
          dscp='BlablaggQWQEQWadweqweqwdqwd123 ajhwlkellkcnl adq'
          socialFB='#!'
          socialLinkedIn='#!'
          socialGithub='#!'
          socialTwitter='#!'
          />
        </Grid>  
        <Grid item xs={cardWidth}>
          <Namecard 
          name="Paella dish" 
          role="CEO, Co-founder" 
          img='https://kiranworkspace.com/demo/projects/code-snippets/team/our-team2/images/users/user1.jpg' 
          dscp='BlablaggQWQEQWadweqweqwdqwd123 ajhwlkellkcnl adq'
          socialFB='#!'
          socialLinkedIn='#!'
          socialGithub='#!'
          socialTwitter='#!'
          />
        </Grid>  
        <Grid item xs={cardWidth}>
          <Namecard 
          name="Paella dish" 
          role="CEO, Co-founder" 
          img='https://kiranworkspace.com/demo/projects/code-snippets/team/our-team2/images/users/user1.jpg' 
          dscp='BlablaggQWQEQWadweqweqwdqwd123 ajhwlkellkcnl adq'
          socialFB='#!'
          socialLinkedIn='#!'
          socialGithub='#!'
          socialTwitter='#!'
          />
        </Grid>  
        <Grid item xs={cardWidth}>
          <Namecard 
          name="Paella dish" 
          role="CEO, Co-founder" 
          img='https://kiranworkspace.com/demo/projects/code-snippets/team/our-team2/images/users/user1.jpg' 
          dscp='BlablaggQWQEQWadweqweqwdqwd123 ajhwlkellkcnl adq'
          socialFB='#!'
          socialLinkedIn='#!'
          socialGithub='#!'
          socialTwitter='#!'
          />
        </Grid>  <Grid item xs={cardWidth}>
          <Namecard 
          name="Paella dish" 
          role="CEO, Co-founder" 
          img='https://kiranworkspace.com/demo/projects/code-snippets/team/our-team2/images/users/user1.jpg' 
          dscp='BlablaggQWQEQWadweqweqwdqwd123 ajhwlkellkcnl adq'
          socialFB='#!'
          socialLinkedIn='#!'
          socialGithub='#!'
          socialTwitter='#!'
          />
        </Grid>  
        
      </Grid>
    )
 }

function Namecard(props) {
  return(
    <Card sx={{ maxWidth: 500,
      display: 'flex',
      paddingLeft:'20px',
      paddingright:'20px'
      }}>
    <CardMedia
      component="img"
      image={props.img}
      alt={props.name}
      sx={{
        borderRadius: '50%',
        height: '150px',
        width: '150px',
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
        <Grid item>
        <IconButton href={props.socialFB}>
          <FacebookIcon/>
        </IconButton>
        </Grid>

        <Grid item>
        <IconButton href={props.socialGithub}>
          <GitHubIcon/>
        </IconButton>
        </Grid>

        <Grid item>
        <IconButton href={props.socialLinkedIn}>
          <LinkedInIcon/>
        </IconButton>
        </Grid>

        <Grid item>
        <IconButton href={props.socialTwitter}>
          <TwitterIcon/>
        </IconButton>
        </Grid>


      </Grid>

    </CardContent>
  </Card>

  );

}
function About() {
  return (
    <div className={styles.headerContainer}>
      <NavBar />
      <Typography sx={{ fontWeight: '400', fontSize: '24px' }}>About us</Typography>
      <text>GeneCruncher was developed by the following team</text>
      
      <MediaCards />


    </div>
  );
}

export default About;