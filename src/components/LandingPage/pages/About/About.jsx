/* eslint-disable jsx-a11y/aria-role */
import React from 'react';
import {
  Typography, Card, CardMedia, CardContent, Grid, IconButton, Divider,
} from '@mui/material';
import FacebookIcon from '@mui/icons-material/Facebook';
import TwitterIcon from '@mui/icons-material/Twitter';
import GitHubIcon from '@mui/icons-material/GitHub';
import LinkedInIcon from '@mui/icons-material/LinkedIn';
import Footer from '../../Footer/Footer';
import NavBar from '../../../NavBar/NavBar';
import styles from './about.module.css';
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

function NameCard(props) {
  const {
    img, name, role, socialFB, socialGithub, socialLinkedIn,
    socialTwitter, dscp,
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
          height: '120px',
          width: '120px',
          margin: '28px',
          objectFit: 'cover',
          border: '4px solid #ff',
          display: 'flex',
        }}
      />

      <CardContent style={{ maxWidth: '200px' }}>
        <Typography variant="h5">{name}</Typography>
        <Typography variant="body2" color="text.secondary">{role}</Typography>
        <br />
        <Typography
          variant="body2"
          style={{
            maxWidth: '200px',
            alignItems: 'center',
            wordWrap: 'break-word',
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
  const { data } = props;
  const cardWidth = 3.5;
  const teamlead = data.shift();
  console.log(teamlead.img);
  return (
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
        sx={{ paddingTop: '50px' }}
      >
        <Grid item xs={cardWidth}>
          <NameCard
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
          <Grid item xs={cardWidth}>
            <NameCard
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
}
);

function About() {
  return (
    <div className={styles.headerContainer}>
      <NavBar />

      <Typography sx={{ fontWeight: '400', fontSize: '24px' }}>Team</Typography>

      {/* For Guy */}
      <Divider variant="middle" textAlign="center" sx={{ padding: '20px' }}>Funder</Divider>

      <Grid
        container
        spacing={0}
        direction="column"
        alignItems="center"
        justifyContent="center"
      >
        <Grid item xs={3}>
          <NameCard
            name="Guy Yachdav"
            role="Supervisor & Initiator"
            img="https://scholar.googleusercontent.com/citations?view_op=medium_photo&user=UoUkGhUAAAAJ&citpid=2"
            dscp="Good looking guy with a lot of experience"
            socialFB=""
            socialLinkedIn="https://www.linkedin.com/in/gyachdav/?originalSubdomain=il"
            socialGithub=""
            socialTwitter=""
          />
        </Grid>

      </Grid>

      <text>GeneCruncher was developed by the following team</text>
      <br />
      <Divider variant="middle" textAlign="left" sx={{ padding: '20px' }}>Frontend Team</Divider>

      <SubteamSection data={frontendData} />

      <Divider variant="middle" textAlign="left" sx={{ padding: '20px' }}>Backend Team</Divider>

      <SubteamSection data={backendData} />

      <Divider variant="middle" textAlign="left" sx={{ padding: '20px' }}>Visualization Team</Divider>

      <SubteamSection data={visualizationData} />

      <Footer />
    </div>
  );
}

export default About;
