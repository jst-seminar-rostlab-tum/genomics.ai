import { Box, IconButton, Typography, Drawer, List, Avatar } from "@mui/material"
import ListIcon from '@mui/icons-material/List';
import { Link, useHistory } from "react-router-dom"

import logo from 'assets/logo.svg';
import { colors } from "shared/theme/colors";
import { useEffect, useRef, useState } from "react";

import { useAuth } from 'shared/context/authContext';
import ProfileImage from "components/ProfileImage";

//In styled(), we cannot use different width to fix different resolution
//we have to use sx
function Appbar(props) {
  return (
    <Box {...props} sx={{
      width: { xs: "90%", sm: "90%", md: "61.8%", lg: "61.8%", xl: "61.8%" },
      display: "flex",
      flexDirection: "row",
      justifyContent: "space-between",
      alignItems: "center",
      backgroundColor: colors.primary[800],
      padding: "0.8em",
      margin: "auto",
      position: "relative"
    }}></Box>
  )
}

function Leftbar(props) {
  return (
    <Box {...props} sx={{
      width: { xs: "100%", sm: "486px", md: "389px", lg: "519px", xl: "664px" },
      display: { xs: "none", sm: "flex", md: "flex", lg: "flex", xl: "flex" },
      flexDirection: "row",
      flexWrap: "warp",
      justifyContent: "space-between",
      alignItems: "center"
    }}
    ></Box>
  )
}

//the DrawerBar will be displayed, if the screen is small, otherwise the Leftbar will be displayed
function DrawerBar({ open, setOpen, executeScroll }) {
  return (
    <Box
      sx={{
        zIndex: "1",
        width: "40%",
        display: { xs: "flex", sm: "none", md: "none", lg: "none", xl: "none" },
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

      <Box sx={{ position: "relative", width: "24px", height: "24px", borderRadius: "5px", bgcolor: colors.primary[400] }}>
        <ListIcon onClick={() => setOpen(!open)} sx={{ position: "absolute", top: "4px", left: "4px", width: "16px", height: "16px", color: "white" }} />

        <Drawer open={open} anchor="bottom" onClose={() => setOpen(false)}>
          <Box sx={{ width: "100vw", height: "25vh", bgcolor: "white", display: "flex", flexDirection: "column", justifyContent: "space-evenly", alignItems: "center" }}>
            <LinkBox to="/explore"><DrawerNavlink>Explore</DrawerNavlink></LinkBox>
            <DrawerNavlink onClick={()=>window.open("https://genecruncher.readthedocs.io/")}>Docs</DrawerNavlink>
            <LinkBox to="/about"><DrawerNavlink>About</DrawerNavlink></LinkBox>
            <Box onClick={executeScroll} sx={{ cursor: "pointer" }}><DrawerNavlink>Contact</DrawerNavlink></Box>
          </Box>
        </Drawer>
      </Box>
    </Box>
  )
}

function Rightbar(props) {
  return (
    <Box {...props} sx={{
      width: { xs: "50%", sm: "216px", md: "166px", lg: "222px", xl: "284px" },
      display: "flex",
      flexDirection: "row",
      alignItems: "center",
      justifyContent: "right",
      gap: "1em"
    }}
    ></Box>
  )
}

function Navlink(props) {

  const { fontWeight="425", main } = props

  return (
    <Box {...props}>
      <Typography {...props} sx={{
        fontWeight,
        fontSize: main ? { xs: "0.6em", sm: "0.8em", md: "0.8em", lg: "1.2em", xl: "1.2em" } : { xs: "0.4em", sm: "0.6em", md: "0.6em", lg: "1em", xl: "1em" },
        color: "white",
        ":hover": {
          color: colors.primary[400],
        }
      }}
      ></Typography>
    </Box>
  )
}

function DrawerNavlink(props) {
  return (
    <Box {...props}>
      <Typography textAlign="center" {...props} sx={{
        fontWeight: "bold",
        fontSize: "1em",
        color: "black",
        ":hover": {
          color: colors.primary[400],
        }
      }}
      ></Typography>
    </Box>
  )
}

function Login(props) {
  return (
    <Box {...props} sx={{
      borderRadius: "10px",
      padding: "8px",
      cursor: "pointer"
    }}
    >
      <Typography {...props} sx={{
        fontWeight: "bold",
        fontSize: { xs: "0.6em", sm: "0.8em", md: "0.8em", lg: "1.2em", xl: "1.2em" },
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

function Signup(props) {
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
        fontSize: { xs: "0.6em", sm: "0.8em", md: "0.8em", lg: "1.2em", xl: "1.2em" },
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

function LinkBox(props) {
  return (
    <Box {...props} component={Link} style={{ textDecoration: "none" }}></Box>
  )
}

export default function Navbar({ 
  onLoginClicked, 
  onSignUpClicked, 
  executeScroll, 
  setNavbarHeight,
  position
}) {

  const [user, setUser] = useAuth()
  const [drawerOpen, setDrawerOpen]=useState(false)
  const history = useHistory()

  function handleClickContactUsInDrawer() {
    console.log("clicked", drawerOpen)
    setDrawerOpen(false)
    console.log("after set", drawerOpen)
    if(position==="fixed") executeScroll()
  }

  //we get the ref of the box that contains the Navbar here
  const boxRef = useRef()

  useEffect(() => {
    //use the set function from Home page to set the height, so that we can use it later
    if(position==="fixed") setNavbarHeight(boxRef.current.clientHeight)
  })

  return (
    <Box ref={boxRef} sx={{width: "100%", bgcolor: colors.primary[800], position: position, zIndex: "10"}}>
      <Appbar>
        <DrawerBar open={drawerOpen} setOpen={setDrawerOpen} executeScroll={handleClickContactUsInDrawer} />
        <Leftbar>
          <LinkBox to="/" sx={{ display: "flex", alignItems: "center", gap: "0.7em" }}>
            <IconButton
              disableRipple
              sx={{
                bgcolor: "white",
                ":hover": { bgcolor: "white" }
              }}
            >
              <img width={28} alt="logo" src={logo} />
            </IconButton>
            <Navlink fontWeight="bold" main={true} >genomics.ai</Navlink>
          </LinkBox>
          <LinkBox to="/explore"><Navlink>Explore</Navlink></LinkBox>
          <Box sx={{cursor: "pointer"}} onClick={()=>window.open("https://genecruncher.readthedocs.io/")}><Navlink>Docs</Navlink></Box>
          <LinkBox to="/about"><Navlink>About</Navlink></LinkBox>
          <Box sx={{ cursor: "pointer" }} onClick={executeScroll}><Navlink>Contact</Navlink></Box>
        </Leftbar>
        <Rightbar>
            {!user && <Login onClick={onLoginClicked}>Log In</Login>}
            {!user && <Signup onClick={onSignUpClicked}>Sign Up</Signup>}
            {
              user && 
              <IconButton onClick={() => history.push("/")}>
                {/* <Avatar alt={user.firstName} src={UserProfileImage} /> */}
                <ProfileImage sizePixels={40} />
              </IconButton>
            }
        </Rightbar>
      </Appbar>
    </Box>
  )
}
