import { Box, IconButton, Typography } from "@mui/material"
import ListIcon from '@mui/icons-material/List';
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
      flexDirection: "row",
      justifyContent: "space-between",
      alignItems: "center",
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
      display: {xs: "none", sm: "flex", md: "flex", lg: "flex", xl: "flex"},
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

function NavlinkDark(props){
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

function Navlink(props){
  return (
    <Box {...props}>
      <Typography {...props} sx={{
        fontWeight: "bold",
        fontSize: {xs: "0.6em", sm: "0.8em", md: "0.8em", lg: "1.2em", xl: "1.2em"},
        color: "black",
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

function DrawerBar({ open, setOpen, executeScroll }){
  return (
      <Box 
        sx={{
          zIndex: "1",
          width: "40%", 
          display: {xs: "flex", sm: "none", md: "none", lg: "none", xl: "none"},
          flexDirection: "row",
          flexWrap: "warp",
          justifyContent: "left",
          alignItems: "center",
          gap: "7%"
        }} 
      >
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
        </LinkBox>

        <Box sx={{position: "relative", width: "24px", height: "24px", borderRadius: "5px", bgcolor: colors.primary[400]}}>
          <ListIcon onClick={()=>setOpen(!open)} sx={{position: "absolute", top: "4px", left: "4px", width: "16px", height: "16px", color: "white"}} />

          <Box 
            sx={{
              position: "absolute", 
              top: "30px", 
              width: "150px", 
              display: open ? "flex" : "none", 
              flexDirection: "column", 
              padding: "5px", 
              bgcolor: "white", 
              borderRadius: "5px", 
              gap: "5px", 
              justifyContent: "space-evenly"
            }}
          >
            <LinkBox to="/about"><Navlink>About us</Navlink></LinkBox>
            <LinkBox to="/docs"><Navlink>Docs </Navlink></LinkBox>
            <Box onClick={executeScroll} sx={{cursor: "pointer"}}><Navlink>Contact us</Navlink></Box>
            <LinkBox to="/explore"><Navlink>Explore</Navlink></LinkBox>
          </Box>
        </Box>
      </Box>
  )
}

export default function Navbar({ onLoginClicked, onSignUpClicked, executeScroll, setNavbarHegiht }) {

  const [drawerOpen, setDrawerOpen]=useState(false)

  function handleClickContactUsInDrawer(){
    console.log("clicked", drawerOpen)
    setDrawerOpen(false)
    console.log("after set", drawerOpen)
    executeScroll()
  }

  const boxRef = useRef()

  useEffect(()=>{
    console.log(boxRef.current.clientHeight)
    setNavbarHegiht(boxRef.current.clientHeight)
  })

  return (
    <Box ref={boxRef} sx={{width: "100%", bgcolor: colors.primary[800], position: "fixed", zIndex: "3"}}>
      <Appbar>
        <DrawerBar open={drawerOpen} setOpen={setDrawerOpen} executeScroll={handleClickContactUsInDrawer} />
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
            <NavlinkDark>genomics.ai</NavlinkDark>
          </LinkBox>
          <LinkBox to="/about"><NavlinkDark>About us</NavlinkDark></LinkBox>
          <LinkBox to="/docs"><NavlinkDark>Docs </NavlinkDark></LinkBox>
          <Box sx={{cursor: "pointer"}} onClick={executeScroll}><NavlinkDark>Contact us</NavlinkDark></Box>
          <LinkBox to="/explore"><NavlinkDark>Explore</NavlinkDark></LinkBox>
        </Leftbar>
        <Rightbar>
          <Login>Log In</Login>
          <Signup>Sign Up</Signup>
        </Rightbar>
      </Appbar>
    </Box>
  )
}