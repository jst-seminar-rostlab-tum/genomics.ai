/* eslint-disable */
import React from 'react';
import {
  Typography, Card, CardMedia, CardContent, Grid, IconButton, Divider, 
} from '@mui/material';
import FacebookIcon from '@mui/icons-material/Facebook';
import TwitterIcon from '@mui/icons-material/Twitter';
import GitHubIcon from '@mui/icons-material/GitHub';
import LinkedInIcon from '@mui/icons-material/LinkedIn';
import Footer from '../../../components/Footer/old';
import NavBar from '../../../components/Navbar/old';
import styles from './about.module.css';
import Stack from "@mui/material/Stack";
import frontendData from './frontend.json';
import backendData from './backend.json';
import visualizationData from './visualization.json';

function FbItem(fblink) {
  return (
    <Grid item>
      <IconButton href={fblink}>
        <FacebookIcon style={{ fill: '#3B8EED' }} />
      </IconButton>
    </Grid>
  );
}

function GithubItem(ghlink) {
  return (
    <Grid item>
      <IconButton href={ghlink}>
        <GitHubIcon style={{ fill: '#171B22' }} />
      </IconButton>
    </Grid>
  );
}

function LinkedInItem(lilink) {
  return (
    <Grid item>
      <IconButton href={lilink}>
        <LinkedInIcon style={{ fill: '#2D69BF' }} />
      </IconButton>
    </Grid>
  );
}

function TwitterItem(ttlink) {
  return (
    <Grid item>
      <IconButton href={ttlink}>
        <TwitterIcon style={{ fill: '#66AFEC' }} />
      </IconButton>
    </Grid>
  );
}

function Namecard(props) {
  const {
    img, name, role, dscp, socialFB, socialGithub, socialLinkedIn,
    socialTwitter,
  } = props;
  return (
    <Card sx={{
      maxWidth: 500,
      display: 'flex',
      paddingLeft: '20px',
      paddingright: '20px',
      boxShadow: 'none',
    }}
    >

      <CardMedia
        component="img"
        image={img}
        alt={name}
        sx={{
          borderRadius: '50%',
          height: '150px',
          width: '150px',
          margin: '28px',
          objectFit: 'cover',
          border: '4px solid #ff',
          display: 'flex',
        }}
      />

      <CardContent style={{ maxWidth: '200px' }}>
        <Typography variant="h5">{name}</Typography>
        <Typography variant="body2" color="text.secondary">{role}</Typography>
        <Typography
          variant="body2"
          style={{
            maxWidth: '200px',
            alignItems: 'center',
            wordWrap: 'break-word',
            paddingTop: '10px'
          }}
        >
          {dscp}
        </Typography>

        <Grid
          container
          spacing={0}
          sx={{
            marginTop: '0px',
            justifyContent: 'center',
          }}
        >

          {
        socialFB !== ''
          ? FbItem(socialFB)
          : <Grid item />
      }
          {
        socialGithub !== ''
          ? GithubItem(socialGithub)
          : <Grid item />
      }
          {
        socialLinkedIn !== ''
          ? LinkedInItem(socialLinkedIn)
          : <Grid item />
      }
          {
        socialTwitter !== ''
          ? TwitterItem(socialTwitter)
          : <Grid item />
      }
        </Grid>

      </CardContent>
    </Card>

  );
}

const SubteamSection = ((props) => {
  const cardWidth = 3.5;
  let { data } = JSON.parse(JSON.stringify(props));
  let teamlead = data.shift();
  return (
    <div >

      <Grid
        container
        spacing={1}
        direction="row"
        justifyContent="center"
        // justifyContent="space-between"
        // alignItems="center"
        // justify="center"
        style={{ minHeight: '20vh', maxWidth: '500vh' }}
        sx={{ paddingTop: '50px' }}
      >
        
          <Grid item xs={cardWidth} style={{ display : 'flex' ,justifyContent: 'center' }}>
            <Namecard
              name={teamlead.name}
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
        sx={{ paddingTop: '50px' }}
      >

        {
        data.map((elem) => (
          <Grid item xs={cardWidth} style={{ display : 'flex' ,justifyContent: 'center' }}>
              <Namecard
                name={elem.name}
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
  );
});

const About = ((props) => {




  return(
    <div>
      <NavBar setUser={props.setUser} />
      <div className={styles.headerContainer}>
        <Stack
          spacing="30px"
        >
          <Typography sx={{ fontWeight: 'bold', fontSize: '30px' }}>Team</Typography>
          <div className={styles.textContainer}>
            <Typography
              sx={{ fontSize: '25px', paddingInline: '350px', maxWidth: '1700px' }}
              align="center"
              >
              Genomics.ai was developed by a team of 12 students from the Technical University of Munich (TUM)
              under the guidance of Dr. Guy Yachdav.
            </Typography>
          </div>

        </Stack>

        <Divider variant="middle" textAlign="left" sx={{ paddingTop: '100px' , paddingBottom: '100px'}}>
          <Typography sx={{ fontSize: '30px', fontWeight: 'bold' }}>Organisation</Typography>
        </Divider>


        <Grid
            container
            spacing={1}
            direction="row"
            justifyContent="center"
            // justifyContent="space-between"
            alignItems="center"
            justify="center"
            style={{ minHeight: '20vh', maxWidth: '500vh' }}
            sx={{ alignItems: 'center' }}
          >
            
            <Grid item xs={3.5} style={{ display : 'flex' ,justifyContent: 'center' }} >
              
                <Namecard
                    name="Dr. Guy Yachdav"
                    role="Supervisor & Initiator"
                    img="https://scholar.googleusercontent.com/citations?view_op=medium_photo&user=UoUkGhUAAAAJ&citpid=2"
                    dscp="Technology executive with over 15 years experience in R&D and specialization in big data and machine learning"
                    socialFB=""
                    socialLinkedIn="https://www.linkedin.com/in/gyachdav/?originalSubdomain=il"
                    socialGithub=""
                    socialTwitter=""
                />
            </Grid>

            
        </Grid>


      <br />
        <Divider variant="middle" textAlign="left" sx={{ paddingTop: '100px' , paddingBottom: '100px'}} >
          <Typography sx={{ fontSize: '30px', fontWeight: 'bold' }}>Frontend</Typography>
        </Divider>

        <SubteamSection data={frontendData} />

        <Divider variant="middle" textAlign="left" sx={{ paddingTop: '100px' , paddingBottom: '100px'}}>
          <Typography sx={{ fontSize: '30px', fontWeight: 'bold' }}>Backend</Typography>
        </Divider>

        <SubteamSection data={backendData} />

        <Divider variant="middle" textAlign="left" sx={{ paddingTop: '100px' , paddingBottom: '100px'}}>
          <Typography sx={{ fontSize: '30px', fontWeight: 'bold' }}>Visualisation</Typography>
        </Divider>
      <SubteamSection data={visualizationData} />

        <Footer />
      </div>
    </div>
  );
});


export default About;
