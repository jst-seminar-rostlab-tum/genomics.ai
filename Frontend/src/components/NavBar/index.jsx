import { Box, IconButton, Typography } from "@mui/material"
import { Link } from "react-router-dom"
import { styled } from "@mui/system";

import logo from 'assets/logo.svg';
import { colors } from "shared/theme/colors";

const Appbar = styled(Box)(({
  width: "100vw",
  display: "flex",
  flexDirection: "row",
  justifyContent: "space-between",
  backgroundColor: colors.primary[800],
  padding: "0.5em",
  paddingLeft: "15%",
  paddingRight: "15%"
})) 

const Leftbar = styled(Box)(({
  width: "70%",
  display: "flex",
  flexDirection: "row",
  justifyContent: "space-between",
  alignItems: "center"
}))

const Rightbar = styled(Box)(({
  width: "30%",
  display: "flex",
  flexDirection: "row",
  alignItems: "center",
  justifyContent: "right",
  gap: "1em"
}))

const Navlink = styled(Typography)(({
  fontWeight: "bold",
  fontSize: "1.2em",
  color: "white",
  ":hover": {
    color: colors.secondary1[200]
  }
}))

const Login = styled(Typography)(({
  fontWeight: "bold",
  fontSize: "1.2em",
  borderRadius: "10px",
  padding: "8px",
  color: "white",
  ":hover": {
    textDecoration: "underline"
  }
}))

const Signup = styled(Typography)(({ theme }) => ({
  borderRadius: "10px",
  border: '2px solid white',
  padding: "8px",
  fontWeight: "bold",
  fontSize: "1.2em",
  color: "white",
  ":hover": {
    border: `2px solid ${theme.palette.primary.light}`,
    color: theme.palette.primary.light,
    textDecoration: "underline"
  }
}))
const Navbar = () => {

  return (
    <Appbar>
      <Leftbar>
        <Link to="/" style={{ textDecoration: "none" }}>
          <Box sx={{ display: "flex", alignItems: "center", gap: "0.7em" }}>
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
          </Box>
        </Link>
        <Link to="/about" style={{ textDecoration: "none" }}>
          <Navlink>
            About us
          </Navlink>
        </Link>
        <Link to="/docs" style={{ textDecoration: "none" }}>
          <Navlink>
           Docs 
          </Navlink>
        </Link>
        <Link to="/contact" style={{ textDecoration: "none" }}>
          <Navlink>
            Contact us
          </Navlink>
        </Link>
        <Link to="/explore" style={{ textDecoration: "none" }}>
          <Navlink>
            Explore
          </Navlink>
        </Link>
      </Leftbar>
      <Rightbar>
        <Link to="/login" style={{ textDecoration: "none", color: "white" }}><Login>Log In</Login></Link>
        <Link to="/signup" style={{ textDecoration: "none", color: "white" }}><Signup>Sign Up</Signup></Link>
      </Rightbar>
    </Appbar>
  )
}

export default Navbar