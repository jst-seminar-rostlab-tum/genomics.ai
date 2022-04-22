import { Box, Typography } from "@mui/material";
import Navbar from "components/Navbar";
import WindowiOS from "components/WindowiOS";
import { useCallback, useState } from "react";
import { colors } from "shared/theme/colors";
import Footer from 'components/Footer';

const Home = ({ setUser }) => {
  const [isLoginFormVisible, setLoginFormVisible] = useState(false);
  const [isRegistrationFormVisible, setRegistrationFormVisible] = useState(false);

  const onLoginClicked = useCallback(() => {
    setRegistrationFormVisible(false)
    setLoginFormVisible(true)
  }, [setLoginFormVisible])

  const onSignUpClicked = useCallback(() => {
    setLoginFormVisible(false);
    setRegistrationFormVisible(true);
  }, [setRegistrationFormVisible])
  
  const onRegistrationFormClosed = useCallback(() => {
    setRegistrationFormVisible(false);
  }, [setRegistrationFormVisible]);

  return (
    <Box sx={{ position: "relative" }}>
      {/* NAVBAR HERE */}
      <Navbar />
      {/* LANDING PAGE */}
      <Box sx={{ 
        position: "absolute",
        backgroundColor: colors.primary[800],
        width: "100vw",
        height: "100vh",
        top: "0",
        left: "0",
        zIndex: "-10"
      }} />
      {/* <Box 
        sx={{
          position: "absolute",
          width: "100vw",
          height: "30vh",
          top: "125%",
          backgroundColor: "white",
          borderRadius: "50%",
          zIndex: "-5"
        }}
      /> */}
      <Box
        sx={{
          width: "100vw",
          paddingLeft: "15%",
          paddingRight: "15%",
        }}
      >
        <WindowiOS />
        <Box>
          <Typography alignText="center" fontSize="2em" fontWeight="bold">What we do</Typography>
        </Box>
      </Box>
    </Box>
  )
}

export default Home