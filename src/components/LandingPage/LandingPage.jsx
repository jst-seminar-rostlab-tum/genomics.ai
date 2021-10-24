import NavBar from '../NavBar/NavBar'
import {useState} from "react";
import LoginForm from "./LoginForm/LoginForm";

function LandingPage() {

    const [isLoginFormVisible, setLoginFormVisible] = useState(false);
    const onLoginClicked = () => {
        setLoginFormVisible(true);
    }
    const onLoginFormClosed = () => {
        setLoginFormVisible(false);
    }

    return (
        <div>
            <NavBar onLoginClicked={onLoginClicked}/>
            <LoginForm visible={isLoginFormVisible} onClose={onLoginFormClosed}/>
        </div>
    );
}

export default LandingPage;