import { Box, IconButton, Typography } from "@mui/material"
import { Link } from "react-router-dom"

import logo from 'assets/logo.svg';
import { colors } from "shared/theme/colors";
import { useEffect, useRef, useState } from "react";

//In styled(), we cannot use different width to fix different resolution
//we have to use sx
function Appbar(props){
  return (
    <Box {...props} sx={{
      width: {xs: "90%", sm: "90%", md: "61.8%", lg: "61.8%", xl: "61.8%"},
      display: "flex",
      flexDirection: {xs: "column", sm: "row", md: "row", lg: "row", xl: "row"},
      justifyContent: "space-between",
      alignItems: {xs: "flex-end", sm: "center", md: "center", lg: "center", xl: "center"},
      backgroundColor: colors.primary[800],
      padding: "0.8em" ,
      margin: "auto",
      position: "relative"
    }}></Box>
  )
}

function Leftbar(props){
  return (
    <Box {...props} sx={{
      width: {xs: "100%", sm: "486px", md: "389px", lg: "519px", xl: "664px"},
      display: "flex",
      flexDirection: "row",
      flexWrap: "warp",
      justifyContent: "space-between",
      alignItems: "center"
    }}
    ></Box>
  )
}

function Rightbar(props){
  return (
    <Box {...props} sx={{
      width: {xs: "50%", sm: "216px", md: "166px", lg: "222px", xl: "284px"},
      display: "flex",
      flexDirection: "row",
      alignItems: "center",
      justifyContent: "right",
      gap: "1em"
    }}
    ></Box>
  )
}

function Navlink(props){
  return (
    <Box {...props}>
      <Typography {...props} sx={{
        fontWeight: "bold",
        fontSize: {xs: "0.6em", sm: "0.8em", md: "0.8em", lg: "1.2em", xl: "1.2em"},
        color: "white",
        ":hover": {
          color: colors.secondary1[200]
        }
      }}
      ></Typography>
    </Box>
  )
}

function Login(props){
  return (
    <Box {...props} sx={{
      borderRadius: "10px",
      padding: "8px",
      cursor: "pointer"
    }}
    >
      <Typography {...props} sx={{
        fontWeight: "bold",
        fontSize: {xs: "0.6em", sm: "0.8em", md: "0.8em", lg: "1.2em", xl: "1.2em"},
        color: "white",
        ":hover": {
          textDecoration: "underline"
        }
      }}
      ></Typography>
    </Box>
  )
}

function Signup(props){
  return (
    <Box {...props} sx={{
      cursor: "pointer",
      borderRadius: "10px",
      border: '2px solid white',
      padding: "8px",
      ":hover": {
        border: `2px solid ${colors.primary[400]}`,
        color: colors.primary[400]
      }
    }}
    >
      <Typography {...props} sx={{
        fontWeight: "bold",
        fontSize: {xs: "0.6em", sm: "0.8em", md: "0.8em", lg: "1.2em", xl: "1.2em"},
        color: "white",
        ":hover": {
          color: colors.primary[400],
          textDecoration: "underline"
        }
      }}
      ></Typography>
    </Box>
  )
}

function LinkBox(props){
  return (
    <Box {...props} component={Link} style={{ textDecoration: "none" }}></Box>
  )
}

export default function Navbar({ onLoginClicked, onSignUpClicked }) {
  return (
    <Appbar>
      <Leftbar>
        <LinkBox to="/" sx={{ display: "flex", alignItems: "center", gap: "0.7em" }}>
          <IconButton
            disableRipple
            sx={{
              bgcolor: "white",
              ":hover": { bgcolor: "primary.dark" }
            }}
            > 
            <img width={28} alt="logo" src={logo} />
          </IconButton>
          <Navlink>genomics.ai</Navlink>
        </LinkBox>
        <LinkBox to="/about"><Navlink>About us</Navlink></LinkBox>
        <LinkBox to="/docs"><Navlink>Docs </Navlink></LinkBox>
        <LinkBox to="/contact"><Navlink>Contact us</Navlink></LinkBox>
        <LinkBox to="/explore"><Navlink>Explore</Navlink></LinkBox>
      </Leftbar>
      <Rightbar>
        <Login>Log In</Login>
        <Signup>Sign Up</Signup>
      </Rightbar>
    </Appbar>
  )
}