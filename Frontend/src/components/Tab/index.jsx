import {useState} from 'react'
import {Link} from 'react-router-dom'

import {Tab, Tabs} from '@mui/material'

/**
 * Styled Tab
 * packaged the Tab from the MUI and add style to fit our mock
 * 
 * need a boolean parameter 
 * - darkBackground
 * in order to render the Tab accordingly
 * 
 * and a label, the content
 * 
 * see the example StyledTab in TabGroup below for how to combine the Tab with a Link from react-router-dom
 */
export function StyledTab(props){

    const {darkBackground}=props

    return(
        <Tab disableRipple {...props}
            sx={{
                textTransform: "none",
                fontWeight: "bold",
                color: darkBackground ? "white" : "black",
                opacity: 1,
                "&.Mui-selected": {
                    fontWeight: "bold",
                    color: darkBackground ? "white" : "black",
                    opacity: 1,
                }
            }}
        />
    )
}

/**
 * Tab Group, the collection of Tabs
 * 
 * like the StyledTab, it also need a boolean parameter 
 * - darkBackground
 * to render the Tab accordingly
 * 
 * it also need an array of objects containing labels and paths
 * 
 * for example: 
 * [
 *     {
 *         label: "tab 1",
 *         path: "/tab_1"
 *     },
 *     {
 *         label: "tab 2",
 *         path: "/tab_2"
 *     },
 *     {
 *         label: "tab 3",
 *         path: "/tab_3"
 *     }
 * ]
 */
export function TabGroup(props) {

    const {darkBackground, tabsInfo}=props
    const [value, setValue]=useState(0)

    return (
        <Tabs value={value} onChange={(_, newValue)=>setValue(newValue)}
            sx={{
                "& .MuiTabs-indicator": {
                    backgroundColor: darkBackground ? "white" : "black"
                }
            }}
        >
            {
                tabsInfo.map((tabInfo)=>(<StyledTab label={tabInfo.label} component={Link} to={tabInfo.path} darkBackground={darkBackground} />))
            }
        </Tabs>
    )
}
