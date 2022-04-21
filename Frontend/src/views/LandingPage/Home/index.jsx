import { Box } from "@mui/material";
import Navbar from "components/Navbar";
import { useCallback, useState } from "react";

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
    <Box
    >
      {/* NAVBAR HERE */}
      <Navbar />
      Hey there
    </Box>
  )
}

export default Home