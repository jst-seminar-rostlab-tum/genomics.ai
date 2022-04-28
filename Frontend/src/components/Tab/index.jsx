import { useRef, useState, useEffect } from 'react'
import { useLocation } from 'react-router-dom'
import { Link } from 'react-router-dom'

import {Tab, Tabs, Box} from '@mui/material'

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

    const position = useLocation()

    const {value, setValue, darkBackground, tabsInfo, width="100%", height="100%"}=props

    const tabsRef = useRef()
    const [tabsHeight, setTabsHeight] = useState(0)

    useEffect(()=>{
        setTabsHeight(tabsRef.current.clientHeight)
    }, [])

    return (
        <Box sx={{width, height, position: "relative"}}>
            <Box ref={tabsRef}>
                <Tabs value={value} onChange={(_, newValue)=>setValue(newValue)}
                    sx={{
                        "& .MuiTabs-indicator": {
                            backgroundColor: darkBackground ? "white" : "black"
                        }
                    }}
                >
                    {
                        tabsInfo.map((tabInfo)=>(<StyledTab label={tabInfo.label} component={Link} to={tabInfo.path ? tabInfo.path : position.pathname} darkBackground={darkBackground} />))
                    }
                </Tabs>
            </Box>
            
            {
                tabsInfo[value].additionalContent && 

                <Box sx={{width: "100%", height: `calc(100% - ${tabsHeight}px)`, overflowY: "scroll"}}>
                    {tabsInfo[value].additionalContent}
                </Box>
            }
        </Box>
    )
}
