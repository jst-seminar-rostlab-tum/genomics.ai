import logo from './logo.svg';
import './app.module.css';
import TopBar from "./components/TopBar/TopBar";
import {Route, BrowserRouter as Router, Switch} from "react-router-dom";
import PageTwo from "./components/PageTwo/PageTwo";
import PageOne from "./components/PageOne/PageOne";
import styles from "./app.module.css"

function App() {
    return (
        <div className={styles.App}>
            <header className={styles.AppHeader}>
                <Router>
                    <TopBar/>
                    <img src={logo} className={styles.AppLogo} alt="logo"/>
                    <Switch>
                        <Route exact path="/" component={PageOne}/>
                        <Route path="/page-two" component={PageTwo}/>
                    </Switch>
                </Router>
            </header>
        </div>
    );
}

export default App;
