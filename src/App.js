import logo from './assets/logo.svg';
import './app.module.css';
import NavigationBar from "./components/NavigationBar/NavigationBar";
import {Route, BrowserRouter as Router, Switch} from "react-router-dom";
import Dashboard from "./components/Dashboard/Dashboard";
import LandingPage from "./components/LandingPage/LandingPage";
import styles from "./app.module.css"

function App() {
    return (
        <div className={styles.App}>
            <header className={styles.AppHeader}>
                <Router>
                    <NavigationBar/>
                    <img src={logo} className={styles.AppLogo} alt="logo"/>
                    <Switch>
                        <Route exact path="/" component={LandingPage}/>
                        <Route path="/page-two" component={Dashboard}/>
                    </Switch>
                </Router>
            </header>
        </div>
    );
}

export default App;
