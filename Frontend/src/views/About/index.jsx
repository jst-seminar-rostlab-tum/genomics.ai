import FacebookIcon from '@mui/icons-material/Facebook';
import TwitterIcon from '@mui/icons-material/Twitter';
import GitHubIcon from '@mui/icons-material/GitHub';
import LinkedInIcon from '@mui/icons-material/LinkedIn';
import { Box, Typography, Grid, IconButton, Divider, Avatar } from '@mui/material';
import { useState, useEffect, useCallback, useRef, useContext } from 'react';
import {
  Switch, Route, Redirect, useHistory, useLocation, useParams,
  useRouteMatch, Link
} from 'react-router-dom';
import NavBar from 'components/NavBar';
import Footer from "components/Footer";
import LoginForm from 'components/LoginForm'
import RegistrationForm from 'components/RegistrationForm'
import { colors } from 'shared/theme/colors';

import organizationData from './data/organization.json'
// import frontendDataOld from './frontend/frontend_old.json';
// import backendDataOld from './backend/backend_old.json';
// import visualizationDataOld from './visualization/visualization_old.json';

import frontendData from './data/frontend.json'
import backendData from './data/backend.json'
import visualizationData from './data/visualization.json'
import { LoginContext } from 'shared/context/loginContext';
import PasswordForgetForm from 'components/PasswordForgetForm';

function FbItem({link}) {
  return (
    <IconButton href={link}>
      <FacebookIcon style={{ fill: '#3B8EED' }} />
    </IconButton>
  );
}

function GithubItem({link}) {
  return (
    <IconButton href={link}>
      <GitHubIcon style={{ fill: '#171B22' }} />
    </IconButton>
  );
}

function LinkedInItem({link}) {
  return (
    <IconButton href={link}>
      <LinkedInIcon style={{ fill: '#2D69BF' }} />
    </IconButton>
  );
}

function TwitterItem({link}) {
  return (
    <IconButton href={link}>
      <TwitterIcon style={{ fill: '#66AFEC' }} />
    </IconButton>
  );
}

function CustomImg({src}){
  return (
    <Avatar alt="Member Avatar" src={src} sx={{width: "80px", height: "80px"}} />
  )
}

function MemberCard(props){

  const fontColor = colors.primary[800]

  const { width="100%", name, roles, img, dscp, socialFB, socialGithub, socialLinkedIn, socialTwitter } = props

  return (
    <Box sx={{width, display: "flex", flexDirection: "column", gap: "0.3em", alignItems: "center"}}>
      <CustomImg src={img} />
      <Typography sx={{color: fontColor, textAlign: "center"}} fontWeight="bold" fontSize="1.2em" >{name}</Typography>
      {
        roles.map((role) => <Typography sx={{color: fontColor, textAlign: "center"}} fontSize="0.9em">{role}</Typography>)
      }
      {/* <Typography sx={{color: fontColor, textAlign: "center"}} fontSize="0.9em">{role}</Typography> */}
      <Typography sx={{color: fontColor, textAlign: "center"}} fontSize="0.9em">{dscp}</Typography>
      <Box sx={{display: "flex", flexDirection: "row", gap: "5px"}}>
        {socialFB ? <FbItem link={socialFB} /> : <></>}
        {socialGithub ? <GithubItem link={socialGithub} /> : <></>}
        {socialLinkedIn ? <LinkedInItem link={socialLinkedIn} /> : <></>}
        {socialTwitter ? <TwitterItem link={socialTwitter} /> : <></>}
      </Box>
    </Box>
  )
}

function MemberSection({name, data}){
  return (
    <Box sx={{marginBottom: "5em", position: "relative", width: { xs: "90%", sm: "90%", md: "70%", lg: "70%", xl: "70%" }}}>
      <Box sx={{ marginBottom: "3em", display: "flex", flexDirection: "row", justifyContent: "space-between", alignItems: "center" }}>
        <Box sx={{ height: "1px", width: { xs: "25%", sm: "37%", md: "35%", lg: "40%", xl: "40%" }, backgroundColor: "black" }} />
        <Typography fontWeight="bold" fontSize={{ xs: "1.2em", sm: "1em", md: "1.4em", lg: "1.7em", xl: "2em" }} >{name}</Typography>
        <Box sx={{ height: "1px", width: { xs: "25%", sm: "37%", md: "35%", lg: "40%", xl: "40%" }, backgroundColor: "black" }} />
      </Box>
      <Grid container rowSpacing={7} columnSpacing={{ xs: 2, sm: 2, md: 3, lg: 4, xl: 4 }}>
        {
          data.map(member => (
            <Grid item xs={12} sm={6} md={4} lg={name==="Organisation" ? 2.4 : 3} xl={name==="Organisation" ? 2.4 : 2} >
              <MemberCard {...member} />
            </Grid>)
          )
        }
      </Grid>
    </Box>
  )
}

export default function About(props){

  const history = useHistory();
  
  const context = useContext(LoginContext)

  const onLoginClicked = () => {
    context.switchRegister(false)
    context.switchLogin(true)
  }

  const onSignUpClicked = () => {
    context.switchLogin(false);
    context.switchRegister(true);
  } 

  const executeScroll = () => history.push({pathname: '/', state: {contact_us: true}})

  return (
    <Box sx={{overflowX: "hidden"}}>
      {context.loginVisible && <LoginForm />}
      {context.registerVisible && <RegistrationForm  />}
      {context.forgetVisible && <PasswordForgetForm />}

      <Box>
        <NavBar
          position="relative"
          onLoginClicked={onLoginClicked}
          onSignUpClicked={onSignUpClicked}
          executeScroll={executeScroll}
        />
      </Box>

      <Box sx={{margin: "3em auto", display: "flex", flexDirection: "column", alignItems: "center", width: { xs: "90%", sm: "90%", md: "70%", lg: "70%", xl: "70%" }}}>
        <Typography fontWeight="bold" fontSize="1.4em">Our Team</Typography>
      </Box>

      <Box sx={{margin: "3em auto", display: "flex", flexDirection: "column", alignItems: "center", width: { xs: "90%", sm: "90%", md: "70%", lg: "70%", xl: "70%" }, textAlign: "center"}}>
        <Typography fontWeight="bold" fontSize="1em">ArchMap was developed by a team of students from the Technical University of Munich (TUM) under the guidance of Dr. Guy Yachdav.</Typography>
      </Box>

      {/* <Box sx={{display: "flex", flexDirection: "column", alignItems: "center"}}>
        <MemberSection name="Organisation" data={organizationData} />
        <MemberSection name="Frontend 1" data={frontend1Data} />
        <MemberSection name="Frontend 2" data={frontend2Data} />
        <MemberSection name="Frontend 3" data={frontend3Data} />
        <MemberSection name="Backend" data={backendData} />
        <MemberSection name="Visualization" data={visualizationData} />
      </Box>

      <Box sx={{margin: "3em auto", display: "flex", flexDirection: "column", alignItems: "center", width: { xs: "90%", sm: "90%", md: "70%", lg: "70%", xl: "70%" }, textAlign: "center"}}>
        <Typography fontWeight="bold" fontSize="1em">Contributors from previous years</Typography>
      </Box> */}

      <Box sx={{display: "flex", flexDirection: "column", alignItems: "center"}}>
        <MemberSection name="Organisation" data={organizationData} />
        <MemberSection name="Frontend" data={frontendData} />
        <MemberSection name="Backend" data={backendData} />
        <MemberSection name="Visualisation" data={visualizationData} />
      </Box>

      <Footer />
    </Box>
  )
}